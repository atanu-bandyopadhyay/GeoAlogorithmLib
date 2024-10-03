#pragma once

#include <math.h>
#include <vector>

//These are dll logic

#define  GEOALGOEXPORT __declspec(dllexport)

using namespace std;
class GEOALGOEXPORT Vector3D
{
public:
	float x;
	float y;
	float z;
	Vector3D() : x(0), y(0), z(0)
	{
	}

	Vector3D(float in_x, float in_y, float in_z) {
		x = in_x;
		y = in_y;
		z = in_z;
	}

	Vector3D(const Vector3D& v3d)
	{
		x = v3d.x;
		y = v3d.y;
		z = v3d.z;
	}

	Vector3D operator-(const Vector3D& nv) const {
		return Vector3D(x - nv.x, y - nv.y, z - nv.z);
	}

	Vector3D cross(const Vector3D& nv) const {
		return Vector3D(
			y * nv.z - z * nv.y,
			z * nv.x - x * nv.z,
			x * nv.y - y * nv.x
		);
	}

	float dot(const Vector3D& nv) const {
		return x * nv.x + y * nv.y + z * nv.z;
	}

	Vector3D normalize() const {
		float len = sqrt(x * x + y * y + z * z);
		return Vector3D(x / len, y / len, z / len);
	}

	~Vector3D()
	{
	}
};

class GEOALGOEXPORT Vector2D
{
public:
	float x;
	float y;
	Vector2D() : x(0), y(0)
	{
	}

	Vector2D(float in_x, float in_y) {
		x = in_x;
		y = in_y;
	}

	Vector2D(const Vector2D& v2d)
	{
		x = v2d.x;
		y = v2d.y;
	}

	~Vector2D()
	{
	}
};

class GEOALGOEXPORT TriSurface {
public:
	Vector3D vt0, vt1, vt2;

	TriSurface(const Vector3D& v0, const Vector3D& v1, const Vector3D& v2)
	{
		vt0 = v0;
		vt1 = v1;
		vt2 = v2;
	}

	TriSurface(const TriSurface& n_tri)
	{
		vt0 = n_tri.vt0;
		vt1 = n_tri.vt1;
		vt2 = n_tri.vt2;
	}

	~TriSurface() {}
};


class GEOALGOEXPORT Line3D
{
public:
	Vector3D startPoint;
	Vector3D endPoint;
	Line3D() {}
	Line3D(Vector3D startP, Vector3D endP)
	{
		startPoint.x = startP.x;
		startPoint.y = startP.y;
		startPoint.z = startP.z;
		endPoint = endP;
		endPoint.x = endP.x;
		endPoint.y = endP.y;
		endPoint.z = endP.z;
	}
	Line3D(const Line3D& l3d)
	{
		startPoint = l3d.startPoint;
		endPoint = l3d.endPoint;
	}
	~Line3D()
	{
	}
};

class GEOALGOEXPORT Line2D
{
public:
	Vector2D startPoint;
	Vector2D endPoint;
	Line2D() {}
	Line2D(Vector2D startP, Vector2D endP)
	{
		startPoint.x = startP.x;
		startPoint.y = startP.y;
		endPoint.x = endP.x;
		endPoint.y = endP.y;
	}
	Line2D(const Line2D& l2d)
	{
		startPoint.x = l2d.startPoint.x;
		startPoint.y = l2d.startPoint.y;
		endPoint.x = l2d.endPoint.x;
		endPoint.y = l2d.endPoint.y;
	}
	~Line2D()
	{
	}
};

float dotProduct(Vector2D p1, Vector2D p2);
GEOALGOEXPORT Vector2D __cdecl getProjectedPointOnLine(Vector2D v1, Vector2D v2, Vector2D p);
GEOALGOEXPORT Vector3D  __cdecl GetVector3D(Vector3D point1, Vector3D point2);
GEOALGOEXPORT Vector2D __cdecl GetVector2D(Vector2D point1, Vector2D point2);
GEOALGOEXPORT Line2D __cdecl parallelLineWithOffset(Line2D lineOrig, float offset, bool opposite = true);
GEOALGOEXPORT bool __cdecl LineLineIntersect2D(Line2D l1, Line2D l2, Vector2D& intersectionPoint);
GEOALGOEXPORT bool __cdecl isPointOnLine(Vector3D p1, Vector3D p2, Vector3D p3);
GEOALGOEXPORT float __cdecl GetDistanceOf2DVector(Vector2D p1, Vector2D p2);
GEOALGOEXPORT float __cdecl GetDistanceOf3DVector(Vector3D p1, Vector3D p2);
GEOALGOEXPORT bool __cdecl isTwoPointSame2D(Vector2D p1, Vector2D p2);
GEOALGOEXPORT float __cdecl GetAngleBetweenTwoLine(Line2D l1, Line2D l2);
std::vector<Line2D> PrepareOffesetPolygon(std::vector<Line2D> polygonDetails, float offset, bool clock = true);
GEOALGOEXPORT bool __cdecl isPointInTriangle(const Vector3D& pt, const TriSurface& tri);
GEOALGOEXPORT Vector3D __cdecl Get3DPointAcoordingToDistance(Vector3D point1, Vector3D normal1, float distance);
GEOALGOEXPORT bool __cdecl getPlaneIntersection(const TriSurface& t1, const TriSurface& t2, Vector3D& intersectionLineStart, Vector3D& intersectionLineEnd);
GEOALGOEXPORT bool __cdecl findIntersectionPoints(const TriSurface& t1, const TriSurface& t2, Vector3D& pt1, Vector3D& pt2);
GEOALGOEXPORT bool __cdecl isTwoPointSame(Vector3D p1, Vector3D p2);
GEOALGOEXPORT bool __cdecl linePlaneIntersection(Vector3D& contact, Vector3D RayPoint2, Vector3D rayOrigin, Vector3D normal, Vector3D planePoint);