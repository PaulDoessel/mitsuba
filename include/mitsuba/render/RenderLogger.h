//
// simple logging api for sampling integrators
//
#pragma once

#include <mitsuba/core/vector.h>
#include <mitsuba/core/normal.h>
#include <mitsuba/core/ray.h>


MTS_NAMESPACE_BEGIN

// basic logger interface
class RenderLogger
{
public:
	virtual ~RenderLogger(){}

	// Logger interface  ------------------------------
	virtual void newPixel( int id )=0;
	virtual void finalizePixel()=0;
	virtual void newSample()=0;
	virtual void newPath(const RayDifferential& ray)=0;
	virtual void pathEvent( const Point3& p )=0;
	virtual void pathEventSurfaceNormal( const Normal& normal )=0;
	virtual void pathEventMesh( const Shape* shape )=0;
	virtual void pathEventLightSample( Spectrum sample, Vector3 direction, float pdf )=0;
	virtual void pathEventBSDFSample( Spectrum sample )=0;
	virtual void pathEventBSDF( Spectrum bsdf, float pdf, int sampledTypeMask )=0;
	virtual void finalizePath()=0;
	virtual void finalizeSample()=0;
};


MTS_NAMESPACE_END