/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/render/RenderLogger.h>
#include <mitsuba/core/statistics.h>

namespace mitsuba
{

static StatsCounter avgPathLength("Path tracer", "Average path length", EAverage);

class PracticalMIPathTracer : public MonteCarloIntegrator
{
public:
	PracticalMIPathTracer(const Properties &props)
		: MonteCarloIntegrator(props)
	{
		m_lightSamples = props.getInteger("lightsamples", 1);
		m_diffuseSamples = props.getInteger("diffusesamples", 1);
		m_glossySamples = props.getInteger("glossysamples", 1);
		m_directLight = props.getBoolean("directlight", true);
		
	}

	/// Unserialize from a binary data stream
	PracticalMIPathTracer(Stream *stream, InstanceManager *manager)
		: MonteCarloIntegrator(stream, manager) { }

	virtual Spectrum Li(const RayDifferential &ray,
		RadianceQueryRecord &rRec ) const
	{
		return Li(ray, rRec, 0);
	}


	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec, RenderLogger* logger) const
	{
		// Some aliases and local variables
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		bool scattered = false;

		if(logger)
			logger->newPath(ray);

		// Perform the first ray intersection (or ignore if the intersection has already been provided).
		rRec.rayIntersect(ray);
		ray.mint = Epsilon;

		Spectrum throughput(1.0f); // this is throughput over pdf
		Float eta = 1.0f;


		if (!its.isValid()) {
			// If no intersection could be found, potentially return
			// radiance from a environment luminaire if it exists
			if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) 
				&& (!m_hideEmitters || scattered))
				Li += throughput * scene->evalEnvironment(ray);
			return Li;
		}
			
		if(logger)
		{
			logger->pathEvent( its.p );
			logger->pathEventSurfaceNormal(its.shFrame.n);
			logger->pathEventMesh( its.shape );
		}

		const BSDF *bsdf = its.getBSDF(ray);

		// Possibly include emitted radiance if requested
		if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
			&& (!m_hideEmitters || scattered))
			Li += throughput * its.Le(-ray.d);


		// Include radiance from a subsurface scattering model if requested
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);


		if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
			|| (m_strictNormals && dot(ray.d, its.geoFrame.n)
				* Frame::cosTheta(its.wi) >= 0)) {

			// Only continue if:
			//   1. The current path length is below the specifed maximum
			//   2. If 'strictNormals'=true, when the geometric and shading
			//      normals classify the incident direction to the same side */
			return Li;
		}


		// ====================================================================
		//                     Direct illumination sampling                    
		// ====================================================================
		// Estimate the direct illumination if this is requested
		DirectSamplingRecord dRec(its);
		if(m_directLight)
		{
		if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
			(bsdf->getType() & BSDF::ESmooth))
		{
			int numSamples = rRec.depth == 1 ? m_lightSamples : 1;

			for(int i=0;i<numSamples;++i)
			{
				Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
				if (!value.isZero())
				{
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					// Allocate a record for querying the BSDF
					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

					// Evaluate BSDF * cos(theta)
					const Spectrum bsdfVal = bsdf->eval(bRec);

					// Prevent light leaks due to the use of shading normals
					if (!bsdfVal.isZero() && (!m_strictNormals
							|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0))
					{
						// Calculate prob. of having generated that direction
						// using BSDF sampling
						Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
							? bsdf->pdf(bRec) : 0;

						// Weight using the power heuristic
						Float weight = miWeight(dRec.pdf, bsdfPdf);
						//Float weight = 1.0f;
						Spectrum lightSample = value * bsdfVal * weight;
						Li += throughput * lightSample/float(numSamples);

						//if(logger)
						//	logger->pathEventLightSample( lightSample, dRec.d, dRec.pdf );
					}
				}
			}
		}

		// ====================================================================
		//                            BSDF sampling                            2222
		// ====================================================================
		{
			// Sample BSDF * cos(theta)
			int numSamples = rRec.depth == 1 ? m_lightSamples : 1;
			for(int i=0;i<numSamples;++i)
			{
				Intersection its_bsdf = its;
				Float bsdfPdf;
				BSDFSamplingRecord bRec(its_bsdf, rRec.sampler, ERadiance);
				Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
				if (bsdfWeight.isZero())
					continue;

				scattered |= bRec.sampledType != BSDF::ENull;

				/// Prevent light leaks due to the use of shading normals
				const Vector wo = its_bsdf.toWorld(bRec.wo);
				Float woDotGeoN = dot(its_bsdf.geoFrame.n, wo);
				if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
					continue;

				bool hitEmitter = false;
				Spectrum value;

				// Trace a ray in this direction
				ray = Ray(its_bsdf.p, wo, ray.time);
				if (scene->rayIntersect(ray, its_bsdf)) {
					// Intersected something - check if it was a luminaire
					if (its_bsdf.isEmitter()) {
						value = its_bsdf.Le(-ray.d);
						dRec.setQuery(ray, its_bsdf);
						hitEmitter = true;
					}
				} else {
					// Intersected nothing -- perhaps there is an environment map?
					const Emitter *env = scene->getEnvironmentEmitter();

					if (env) {
						value = env->evalEnvironment(ray);
						if (!env->fillDirectSamplingRecord(dRec, ray))
							continue;
						hitEmitter = true;
					} else {
						continue;
					}
				}

				// If a luminaire was hit, estimate the local illumination and
				// weight using the power heuristic
				if (hitEmitter &&
					(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance))
				{
					// Compute the prob. of generating that direction using the
					// implemented direct illumination sampling technique
					const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
						scene->pdfEmitterDirect(dRec) : 0;
					Spectrum bsdfSample = throughput * bsdfWeight * value * miWeight(bsdfPdf, lumPdf);
					Li += bsdfSample/float(numSamples);
				}
			}
		} // bsdf samples
		} // direct light

		// ====================================================================
		//                         Indirect illumination                       
		// ====================================================================

		///*
		// indirect diffuse -------
		if(its.shape->getBSDF()->getType() & BSDF::EDiffuse)
		{
			// Sample BSDF * cos(theta)
			int numSamples = rRec.depth == 1 ? m_diffuseSamples : 1;
			RadianceQueryRecord indirectDiffuseRec;
			for(int i=0;i<numSamples;++i)
			{
				Intersection diffuse_its = its;
				//std::cout << "indirect diffuse sample\n"; std::flush(std::cout);
				Float bsdfPdf;
				BSDFSamplingRecord bRec(diffuse_its, rRec.sampler, ERadiance);
				bRec.typeMask = BSDF::EDiffuse;
				Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
				if (bsdfWeight.isZero())
					continue;

				indirectDiffuseRec.recursiveQuery(rRec, RadianceQueryRecord::ERadianceNoEmission);

				const Vector wo = diffuse_its.toWorld(bRec.wo);
				RayDifferential indirectDiffuseRay(diffuse_its.p, wo, ray.time);

				Spectrum indirectDiffuseSample = throughput * bsdfWeight * Li_path(indirectDiffuseRay, indirectDiffuseRec, 0);
				Li += indirectDiffuseSample/float(numSamples);
			}
		} // indirect diffuse
		// indirect glossy -------
		if(its.shape->getBSDF()->getType() & BSDF::EGlossy)
		{
			// Sample BSDF * cos(theta)
			int numSamples = rRec.depth == 1 ? m_glossySamples : 1;
			RadianceQueryRecord indirectGlossyRec;
			for(int i=0;i<numSamples;++i)
			{
				Intersection glossy_its = its;
				//std::cout << "indirect diffuse sample\n"; std::flush(std::cout);
				Float bsdfPdf;
				BSDFSamplingRecord bRec(glossy_its, rRec.sampler, ERadiance);
				bRec.typeMask = BSDF::EGlossy;
				Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
				if (bsdfWeight.isZero())
					continue;

				indirectGlossyRec.recursiveQuery(rRec, RadianceQueryRecord::ERadianceNoEmission);

				const Vector wo = glossy_its.toWorld(bRec.wo);
				RayDifferential indirectDiffuseRay(glossy_its.p, wo, ray.time);

				Spectrum indirectDiffuseSample = throughput * bsdfWeight * Li_path(indirectDiffuseRay, indirectGlossyRec, 0);
				Li += indirectDiffuseSample/float(numSamples);
			}
		} // indirect diffuse
		//*/

		if(logger)
			logger->finalizePath();

		// Store statistics
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;

		return Li;
	}

	Spectrum Li_path(const RayDifferential &r, RadianceQueryRecord &rRec, RenderLogger* logger) const
	{
		// Some aliases and local variables
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		bool scattered = true;


		// Perform the first ray intersection (or ignore if the intersection has already been provided).
		rRec.rayIntersect(ray);
		ray.mint = Epsilon;

		Spectrum throughput(1.0f); // this is throughput over pdf
		Float eta = 1.0f;

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0)
		{
			//std::cout << "indirect diffuse sample-check0\n"; std::flush(std::cout);
			if (!its.isValid()) {
				// If no intersection could be found, potentially return
				// radiance from a environment luminaire if it exists
				if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) 
					&& (!m_hideEmitters || scattered))
					Li += throughput * scene->evalEnvironment(ray);
				break;
			}
			
			if(logger)
			{
				logger->pathEvent( its.p );
				logger->pathEventSurfaceNormal(its.shFrame.n);
				logger->pathEventMesh( its.shape );
			}

			const BSDF *bsdf = its.getBSDF(ray);

			// Possibly include emitted radiance if requested
			if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
				&& (!m_hideEmitters || scattered))
				Li += throughput * its.Le(-ray.d);


			// Include radiance from a subsurface scattering model if requested
			if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
				Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);


			if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
				|| (m_strictNormals && dot(ray.d, its.geoFrame.n)
					* Frame::cosTheta(its.wi) >= 0)) {

				// Only continue if:
				//   1. The current path length is below the specifed maximum
				//   2. If 'strictNormals'=true, when the geometric and shading
				//      normals classify the incident direction to the same side */
				break;
			}


			// ====================================================================
			//                     Direct illumination sampling                    
			// ====================================================================
			// Estimate the direct illumination if this is requested
			DirectSamplingRecord dRec(its);
			///*
			if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
				(bsdf->getType() & BSDF::ESmooth))
			{
				Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
				if (!value.isZero())
				{
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					// Allocate a record for querying the BSDF
					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

					// Evaluate BSDF * cos(theta)
					const Spectrum bsdfVal = bsdf->eval(bRec);

					// Prevent light leaks due to the use of shading normals
					if (!bsdfVal.isZero() && (!m_strictNormals
							|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0))
					{
						// Calculate prob. of having generated that direction
						// using BSDF sampling
						Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
							? bsdf->pdf(bRec) : 0;

						// Weight using the power heuristic
						Float weight = miWeight(dRec.pdf, bsdfPdf);
						Spectrum lightSample = value * bsdfVal * weight;
						Li += throughput * lightSample;

						//if(logger)
						//	logger->pathEventLightSample( lightSample, dRec.d, dRec.pdf );
					}
				}
			}
			//*/

			///*
			// ====================================================================
			//                            BSDF sampling                            
			// ====================================================================

			// Sample BSDF * cos(theta)
			Float bsdfPdf;
			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
			Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
			if (bsdfWeight.isZero())
				break;

			if(logger)
				logger->pathEventBSDF( bsdfWeight*bsdfPdf, bsdfPdf, bRec.sampledType );

			scattered |= bRec.sampledType != BSDF::ENull;

			/// Prevent light leaks due to the use of shading normals
			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				break;

			bool hitEmitter = false;
			Spectrum value;

			// Trace a ray in this direction
			ray = Ray(its.p, wo, ray.time);
			if (scene->rayIntersect(ray, its)) {
				// Intersected something - check if it was a luminaire
				if (its.isEmitter()) {
					value = its.Le(-ray.d);
					dRec.setQuery(ray, its);
					hitEmitter = true;
				}
			} else {
				// Intersected nothing -- perhaps there is an environment map?
				const Emitter *env = scene->getEnvironmentEmitter();

				if (env) {
					if (m_hideEmitters && !scattered)
						break;

					value = env->evalEnvironment(ray);
					if (!env->fillDirectSamplingRecord(dRec, ray))
						break;
					hitEmitter = true;
				} else {
					break;
				}
			}

			// Keep track of the throughput and relative
			// refractive index along the path
			throughput *= bsdfWeight;
			eta *= bRec.eta;

			// If a luminaire was hit, estimate the local illumination and
			// weight using the power heuristic
			if (hitEmitter &&
				(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
				// Compute the prob. of generating that direction using the
				// implemented direct illumination sampling technique
				const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
					scene->pdfEmitterDirect(dRec) : 0;
				Spectrum bsdfSample = throughput * value * miWeight(bsdfPdf, lumPdf);
				Li += bsdfSample;

				if(logger)
					logger->pathEventBSDFSample( bsdfSample );
			}

			// ====================================================================
			//                         Indirect illumination                       
			// ====================================================================

			// Set the recursive query type. Stop if no surface was hit by the
			// BSDF sample or if indirect illumination was not requested
			if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
				break;
			rRec.type = RadianceQueryRecord::ERadianceNoEmission;

			if (rRec.depth++ >= m_rrDepth) {
				// Russian roulette: try to keep path weights equal to one,
				// while accounting for the solid angle compression at refractive
				// index boundaries. Stop with at least some probability to avoid
				// getting stuck (e.g. due to total internal reflection)
				Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
				if (rRec.nextSample1D() >= q)
					break;
				throughput /= q;
			}
			//*/
		}

		// Store statistics
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;

		return Li;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MIPathTracer[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}


	int m_lightSamples;
	int m_diffuseSamples;
	int m_glossySamples;
	bool m_directLight;

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(PracticalMIPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(PracticalMIPathTracer, "practical MI path tracer");
} // namespace mitsuba

