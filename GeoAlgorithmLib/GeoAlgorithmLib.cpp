/* 
 * GeoAlgorithmLib - Indranil Ghosh Roy, Atanu Bandyopadhyay, Subhadeep Kanungo
 * 
 * Copyright (C) 2023-2024 EPFL SCI STI MM
 *
 * This file is part of GeoAlgorithmLib.
 *
 * GeoAlgorithmLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GeoAlgorithmLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GeoAlgorithmLib.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Additional permission under GNU GPL version 3 section 7
 * 
 * If you modify this Program, or any covered work, by linking or combining it
 * with Eclipse (or a modified version of Eclipse or an Eclipse plugin or 
 * an Eclipse library), containing parts covered by the terms of the 
 * Eclipse Public License (EPL), the licensors of this Program grant you 
 * additional permission to convey the resulting work.  Corresponding Source 
 * for a non-source form of such a combination shall include the source code 
 * for the parts of Eclipse libraries used as well as that of the  covered work.
 * 
 */
#include "pch.h"
#include "GeoAlgorithmLib.hpp"

#define epsilon_intersection 0.001

float dotProduct(Vector2D p1, Vector2D p2)
{
    return (p1.x * p2.x + p1.y * p2.y);
}

Vector2D getProjectedPointOnLine(Vector2D v1, Vector2D v2, Vector2D p)
{
    // get dot product of e1, e2
    Vector2D e1(v2.x - v1.x, v2.y - v1.y);
    Vector2D e2(p.x - v1.x, p.y - v1.y);
    double valDp = dotProduct(e1, e2);
    // get length of vectors
    double lenLineE1 = sqrt(e1.x * e1.x + e1.y * e1.y);
    double lenLineE2 = sqrt(e2.x * e2.x + e2.y * e2.y);
    double cos = valDp / (lenLineE1 * lenLineE2);
    // length of v1P'
    double projLenOfLine = cos * lenLineE2;
    Vector2D p1((v1.x + (projLenOfLine * e1.x) / lenLineE1),
        (v1.y + (projLenOfLine * e1.y) / lenLineE1));
    return p1;
}

Vector3D GetVector3D(Vector3D point1, Vector3D point2)
{
    Vector3D vector;
    vector.x = point1.x - point2.x;
    vector.y = point1.y - point2.y;
    vector.z = point1.z - point2.z;
    return vector;
}

Vector2D GetVector2D(Vector2D point1, Vector2D point2)
{
    Vector2D vector;
    vector.x = point1.x - point2.x;
    vector.y = point1.y - point2.y;
    return vector;
}

Line2D parallelLineWithOffset(Line2D lineOrig, float offset, bool opposite)
{
    float dx = lineOrig.endPoint.x - lineOrig.startPoint.x;
    float dy = lineOrig.endPoint.y - lineOrig.startPoint.y;

    //a vector perpendicular to EndPoint - StartPoint
    float peX = dy;
    float peY = -dx;

    if (opposite) {
        peX = -dy;
        peY = dx;
    }

    float len = sqrt(peX * peX + peY * peY);
    peX /= len;
    peY /= len;

    peX *= offset;
    peY *= offset;

    Line2D parallelLine;

    parallelLine.startPoint.x = lineOrig.startPoint.x + peX;
    parallelLine.startPoint.y = lineOrig.startPoint.y + peY;

    parallelLine.endPoint.x = lineOrig.endPoint.x + peX;
    parallelLine.endPoint.y = lineOrig.endPoint.y + peY;
    return parallelLine;
}

GEOALGOEXPORT bool isTwoPointSame2D(Vector2D p1, Vector2D p2) {
    if ((fabs(p1.x - p2.x) < epsilon_intersection) && (fabs(p1.y - p2.y) < epsilon_intersection)) {
        return true;
    }
    return false;
}

float GetDistanceOf2DVector(Vector2D p1, Vector2D p2)
{
    float distance = 0.0;
    float x0 = p2.x;
    float y0 = p2.y;
    float x1 = p1.x;
    float y1 = p1.y;
    distance = (float)(sqrt((double)(((x0 - x1) * (x0 - x1)) + ((y0 - y1) * (y0 - y1)))));
    return distance;
}


float Det(float a, float b, float c, float d)
{
    return a * d - b * c;
}

float GetAngleBetweenTwoLine(Line2D l1, Line2D l2) {
    float a = l1.endPoint.x - l1.startPoint.x;
    float b = l1.endPoint.y - l1.startPoint.y;
    //float c = l2.startPoint.x - l2.endPoint.x;
    //float d = l2.startPoint.y - l2.endPoint.y;
    float c = l2.endPoint.x - l2.startPoint.x;
    float d = l2.endPoint.y - l2.startPoint.y;

    //
    float cos_angle, angle;
    float mag_v1 = sqrt(a * a + b * b);
    float mag_v2 = sqrt(c * c + d * d);
    //
    cos_angle = (a * c + b * d) / (mag_v1 * mag_v2);
    cos_angle = fmax(-1.0, fmin(1.0, cos_angle));
    angle = acos(cos_angle);
    angle = angle * (180.0 / 3.14159); // convert to degrees
    //

    return angle;
}

bool LineLineIntersect2D(Line2D l1, Line2D l2, Vector2D& intersectionPoint)
{
    if (isTwoPointSame2D(l1.endPoint, l2.startPoint))
    {
        intersectionPoint = l1.endPoint;
        return true;
    }

    float d1 = GetDistanceOf2DVector(l2.startPoint, l2.endPoint);
    float d2 = GetDistanceOf2DVector(l1.endPoint, l2.startPoint);
    float d3 = GetDistanceOf2DVector(l1.endPoint, l2.endPoint);

    if (fabs(d1 - (d2 + d3)) < 0.01) {
        intersectionPoint = l1.endPoint;
        return true;
    }


    float x1 = l1.startPoint.x;
    float y1 = l1.startPoint.y;

    float x2 = l1.endPoint.x;
    float y2 = l1.endPoint.y;

    float x3 = l2.startPoint.x;
    float y3 = l2.startPoint.y;

    float x4 = l2.endPoint.x;
    float y4 = l2.endPoint.y;


    float detL1 = Det(x1, y1, x2, y2);
    float detL2 = Det(x3, y3, x4, y4);
    float x1mx2 = x1 - x2;
    float x3mx4 = x3 - x4;
    float y1my2 = y1 - y2;
    float y3my4 = y3 - y4;

    float xnom = Det(detL1, x1mx2, detL2, x3mx4);
    float ynom = Det(detL1, y1my2, detL2, y3my4);
    float denom = Det(x1mx2, y1my2, x3mx4, y3my4);
    //Lines are parallel
    if (denom == 0.0)
    {
        intersectionPoint.x = NAN;
        intersectionPoint.y = NAN;
        return false;
    }

    intersectionPoint.x = xnom / denom;
    intersectionPoint.y = ynom / denom;
    if (!isfinite(intersectionPoint.x) || !isfinite(intersectionPoint.y)) //Probably a numerical issue
        return false;

    return true;
}

std::vector<Line2D> PrepareOffesetPolygon(std::vector<Line2D> polygonDetails, float offset, bool clock)
{
    std::vector<Line2D> results;
    if (polygonDetails.size() < 3)
    {
        return polygonDetails;
    }

    for (int i = 0; i < polygonDetails.size(); i++)
    {
        Line2D line2D = polygonDetails[i];
        float lineLength = GetDistanceOf2DVector(line2D.startPoint, line2D.endPoint);
        if (lineLength < offset) {
            if (i < polygonDetails.size() - 1)
            {
                polygonDetails[i + 1].startPoint = polygonDetails[i].startPoint;
            }
            continue;
        }
        if (i < polygonDetails.size() - 1)
        {
            Line2D line2DNext = polygonDetails[i + 1];
            float angle = GetAngleBetweenTwoLine(line2D, line2DNext);
            if (fabs(angle - 180) < 1 || fabs(angle) < 1 || fabs(angle - 360) < 1)
            {
                polygonDetails[i + 1].startPoint = polygonDetails[i].startPoint;
                continue;
            }

        }

        Line2D parallelLine2D = parallelLineWithOffset(line2D, offset, clock);
        results.push_back(parallelLine2D);

        if (results.size() > 1)
        {
            Line2D l1 = results[results.size() - 2];
            Line2D l2 = results[results.size() - 1];
            Vector2D intersectionPoint;
            bool val1 = LineLineIntersect2D(l1, l2, intersectionPoint);
            if (val1) {
                float dist1 = GetDistanceOf2DVector(results[results.size() - 2].endPoint, intersectionPoint);
                float dist2 = GetDistanceOf2DVector(results[results.size() - 1].startPoint, intersectionPoint);
                if (dist1 > (offset * 2) || (dist2 > offset * 2)) {
                    results[results.size() - 2].endPoint = results[results.size() - 1].startPoint;
                }
                else {
                    results[results.size() - 2].endPoint = intersectionPoint;
                    results[results.size() - 1].startPoint = intersectionPoint;
                }
            }
        }
    }
    if (results.size() > 0) {
        results.erase(results.begin());
    }
    return results;
}

GEOALGOEXPORT float GetDistanceOf3DVector(Vector3D p1, Vector3D p2)
{
    float distance = 0.0;
    float x0 = p2.x;
    float y0 = p2.y;
    float z0 = p2.z;
    float x1 = p1.x;
    float y1 = p1.y;
    float z1 = p1.z;
    distance = (float)(sqrt((double)(((x0 - x1) * (x0 - x1)) + ((y0 - y1) * (y0 - y1)) + ((z0 - z1) * (z0 - z1)))));
    return distance;
}

GEOALGOEXPORT bool isPointOnLine(Vector3D p1, Vector3D p2, Vector3D p3) {
    float d1 = GetDistanceOf3DVector(p1, p2);
    float d2 = GetDistanceOf3DVector(p1, p3);
    float d3 = GetDistanceOf3DVector(p2, p3);
    if (fabs(d3 - (d1 + d2)) < 0.1) {
        return true;
    }
    return false;
}

float triangleArea(const Vector3D& a, const Vector3D& b, const Vector3D& c) {
    Vector3D ab = b - a;
    Vector3D ac = c - a;
    Vector3D norm = ab.cross(ac).normalize();
    return 0.5f * norm.x + 0.5f * norm.y + 0.5f * norm.y; // Area = 0.5 * |AB x AC|
}

bool isPointInTriangleArea(const Vector3D& p, const Vector3D& a, const Vector3D& b, const Vector3D& c) {
    float areaABC = triangleArea(a, b, c);
    float areaABP = triangleArea(a, b, p);
    float areaACP = triangleArea(a, c, p);
    float areaBCP = triangleArea(b, c, p);

    // The point is inside the triangle if the sum of areas of the sub-triangles equals the area of the triangle
    return (areaABP + areaACP + areaBCP) <= areaABC;
}

Vector3D calculateNormal(const Vector3D& a, const Vector3D& b, const Vector3D& c) {
    Vector3D ab = b - a;
    Vector3D ac = c - a;
    return ab.cross(ac).normalize(); // Normal vector of the triangle
}

bool isSameSide(const Vector3D& p, const Vector3D& a, const Vector3D& b, const Vector3D& c) {
    Vector3D normal = calculateNormal(a, b, c);
    float d = normal.dot(a);
    float planeValueP = normal.dot(p) - d;

    if (planeValueP > 0) return false;  // P is on the positive side of the plane
    // Check if point P is on the same side of each triangle edge


    Vector3D ab = b - a;
    Vector3D ap = p - a;
    Vector3D ac = c - a;

    // Using cross products to determine the side of the edges
    Vector3D cross1 = normal.cross(ab); // Cross product of normal and edge AB
    Vector3D cross2 = normal.cross(ac); // Cross product of normal and edge AC
    Vector3D cross3 = normal.cross(c - b); // Cross product of normal and edge BC

    float side1 = cross1.dot(ap);
    float side2 = cross2.dot(p - b);
    float side3 = cross3.dot(p - c);

    // Point is inside if all dot products are non-negative (or non-positive)
    return ((side1 >= 0 && side2 >= 0 && side3 >= 0) || (side1 <= 0 && side2 <= 0 && side3 <= 0));
}

bool isPointInTriangleNew(const Vector3D& p, const Vector3D& a, const Vector3D& b, const Vector3D& c) {
    Vector3D v0 = b - a;
    Vector3D v1 = c - a;
    Vector3D v2 = p - a;

    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);

    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    // Check if the point is in the triangle
    return (u >= 0) && (v >= 0) && (w >= 0);
}


GEOALGOEXPORT bool isPointInTriangle(const Vector3D& pt, const TriSurface& tri) {
    //return isPointInTriangleNew(pt, tri.vt0, tri.vt1, tri.vt2);
    //return isSameSide(pt, tri.vt0, tri.vt1, tri.vt2);
    Vector3D v0 = tri.vt1 - tri.vt0;
    Vector3D v1 = tri.vt2 - tri.vt0;
    Vector3D v2 = pt - tri.vt0;

    float dot00 = v0.dot(v0);
    float dot01 = v0.dot(v1);
    float dot02 = v0.dot(v2);
    float dot11 = v1.dot(v1);
    float dot12 = v1.dot(v2);

    float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u >= 0) && (v >= 0) && (u + v <= 1);
}

GEOALGOEXPORT Vector3D Get3DPointAcoordingToDistance(Vector3D point1, Vector3D normal1, float distance)
{
    Vector3D newPoint3D;
    newPoint3D.x = point1.x + distance * normal1.x;
    newPoint3D.y = point1.y + distance * normal1.y;
    newPoint3D.z = point1.z + distance * normal1.z;
    return newPoint3D;
}

GEOALGOEXPORT bool isTwoPointSame(Vector3D p1, Vector3D p2) {
    if (fabs(p1.x - p2.x) < 0.001 && fabs(p1.y - p2.y) < 0.001 && fabs(p1.z - p2.z) < 0.001) {
        return true;
    }
    return false;
}

GEOALGOEXPORT bool getPlaneIntersection(const TriSurface& t1, const TriSurface& t2, Vector3D& intersectionLineStart, Vector3D& intersectionLineEnd) {
    Vector3D e1 = t1.vt1 - t1.vt0;
    Vector3D e2 = t1.vt2 - t1.vt0;
    Vector3D n1 = e1.cross(e2);

    Vector3D e3 = t2.vt1 - t2.vt0;
    Vector3D e4 = t2.vt2 - t2.vt0;
    Vector3D n2 = e3.cross(e4);

    Vector3D d = t2.vt0 - t1.vt0;
    Vector3D r = n1.cross(n2);

    float denom = r.dot(r);
    if (fabs(denom) < 1e-6) {
        return false; // Planes are parallel
    }

    float t = d.dot(n2.cross(n1)) / denom;
    Vector3D origin;
    origin.x = t1.vt0.x;
    origin.x = t1.vt0.y;
    origin.x = t1.vt0.z;
    intersectionLineStart = Get3DPointAcoordingToDistance(origin, e1, t);
    intersectionLineEnd = Get3DPointAcoordingToDistance(origin, e2, t);

    return true;
}

GEOALGOEXPORT Vector3D getPlaneNormalFromThreePoints(Vector3D p1, Vector3D p2, Vector3D p3)
{
    Vector3D dirX = p2 - p1;
    Vector3D dirY = p3 - p1;
    Vector3D dirZ = dirX.cross(dirY);
    Vector3D unitZ = dirZ.normalize();
   // unitZ.x = -unitZ.x;
   // unitZ.y = -unitZ.y;
   // unitZ.z = -unitZ.z;
    return unitZ;
}

GEOALGOEXPORT bool linePlaneIntersection(Vector3D& contact, Vector3D RayPoint2, Vector3D rayOrigin, Vector3D normal, Vector3D planePoint) {

    Vector3D v = RayPoint2 - rayOrigin;
    Vector3D w = planePoint - rayOrigin;

    float d = normal.dot(v);
    //if ( d == 0){
    //	return false;
    //}

    float k = normal.dot(w) / d;

    contact.x = rayOrigin.x + (k * v.x);//normalV.x);
    contact.y = rayOrigin.y + (k * v.y);//(k*normalV.y);
    contact.z = rayOrigin.z + (k * v.z);//(k*normalV.z);

   // if (fabs(planePoint.z - contact.z) < 0.1) {
        Vector3D cp;
        cp.x = contact.x;
        cp.y = contact.y;
        cp.z = contact.z;

        float d1 = GetDistanceOf3DVector(rayOrigin, RayPoint2);
        float d2 = GetDistanceOf3DVector(rayOrigin, cp);
        float d3 = GetDistanceOf3DVector(cp, RayPoint2);

        if (fabs(d1 - (d2 + d3)) < 0.001) {
            return true;
        }
    //}
    return false;
}

GEOALGOEXPORT bool getPlaneIntersectionPoints(Vector3D planeOrigin, Vector3D planeNormal, TriSurface t1, Vector3D& l1_start, Vector3D& l1_end) {
    Vector3D tp11;
    tp11.x = t1.vt0.x;
    tp11.y = t1.vt0.y;
    tp11.z = t1.vt0.z;

    Vector3D tp12;
    tp12.x = t1.vt1.x;
    tp12.y = t1.vt1.y;
    tp12.z = t1.vt1.z;

    Vector3D tp13;
    tp13.x = t1.vt2.x;
    tp13.y = t1.vt2.y;
    tp13.z = t1.vt2.z;

    Vector3D intersectionP1;
    Vector3D intersectionP2;
    Vector3D intersectionP3;
    bool f1 = linePlaneIntersection(intersectionP1, tp11, tp12, planeNormal, planeOrigin);
    bool f2 = linePlaneIntersection(intersectionP2, tp11, tp13, planeNormal, planeOrigin);
    bool f3 = linePlaneIntersection(intersectionP3, tp13, tp12, planeNormal, planeOrigin);



    if ((!f1 && !f2) || (!f2 && !f3) || (!f1 && !f3))
    {
        return false;
    }

    if (f1 && f2)
    {
        if (!isTwoPointSame(intersectionP1, intersectionP2)) {
            l1_start = intersectionP1;
            l1_end = intersectionP2;
        }
        //else {
        //    f2 = false;
        //}
    }

    if (f1 && f3)
    {
        if (!isTwoPointSame(intersectionP1, intersectionP3)) {
            l1_start = intersectionP1;
            l1_end = intersectionP3;
        }
        //else {
        //    f3 = false;
        //}
    }

    if (f2 && f3)
    {
        if (!isTwoPointSame(intersectionP2, intersectionP3)) {
            l1_start = intersectionP2;
            l1_end = intersectionP3;
        }
        //else {
        //    f3 = false;
        //}
    }
    return true;
}

bool isPointOnEdge(Vector3D pt, TriSurface t1) {
    //Line3D l1
    bool f1 = isPointOnLine(pt, t1.vt0, t1.vt1);
    bool f2 = isPointOnLine(pt, t1.vt1, t1.vt2);
    bool f3 = isPointOnLine(pt, t1.vt2, t1.vt0);

    if (f1 || f2 || f3)
    {
        return true;
    }
    return false;
}

GEOALGOEXPORT bool findIntersectionPoints(const TriSurface& t1, const TriSurface& t2, Vector3D& pt1, Vector3D& pt2) {
    bool p1_found = false;
    bool p2_found = false;
    Vector3D l1_start, l1_end;

    Vector3D tp21;
    tp21.x = t2.vt0.x;
    tp21.y = t2.vt0.y;
    tp21.z = t2.vt0.z;

    Vector3D tp22;
    tp22.x = t2.vt1.x;
    tp22.y = t2.vt1.y;
    tp22.z = t2.vt1.z;

    Vector3D tp23;
    tp23.x = t2.vt2.x;
    tp23.y = t2.vt2.y;
    tp23.z = t2.vt2.z;

    Vector3D tp11;
    tp11.x = t1.vt0.x;
    tp11.y = t1.vt0.y;
    tp11.z = t1.vt0.z;

    Vector3D tp12;
    tp12.x = t1.vt1.x;
    tp12.y = t1.vt1.y;
    tp12.z = t1.vt1.z;

    Vector3D tp13;
    tp13.x = t1.vt2.x;
    tp13.y = t1.vt2.y;
    tp13.z = t1.vt2.z;

    Vector3D normalV = getPlaneNormalFromThreePoints(tp21, tp22, tp23);
    bool firstTime = true;
    bool secondTime = true;

    if (!getPlaneIntersectionPoints(tp21, normalV, t1, l1_start, l1_end))
    {
        firstTime = false;
        //return false;
    }

    //if (isTwoPointSame(l1_start, l1_end)) {
    //    return false;
    //}


    // Check intersection points with the edges of the second triangle
    // TODO: Implement edge intersection checks
    //if (!isPointInTriangle(l1_start, t2) || !isPointInTriangle(l1_end, t2))
    //{
        //recalculate with 2nd triangle 
    Vector3D l2_start, l2_end;
    Vector3D normalU = getPlaneNormalFromThreePoints(tp11, tp12, tp13);

    
    if (!getPlaneIntersectionPoints(tp11, normalU, t2, l2_start, l2_end))
    {
        secondTime = false;
        //return false;
    }

    if (!firstTime && !secondTime)
    {
        return false;
    }
        
    bool startFound = false;
    bool endFound = false;
    Vector3D newStart;
    Vector3D newEnd;

    if (firstTime) {
        if (isPointInTriangle(l1_start, t1) && isPointInTriangle(l1_start, t2))
        {
            startFound = true;
            newStart = l1_start;
        }

        if (!startFound) {
            if (isPointInTriangle(l1_start, t1))
            {
                if (isPointOnEdge(l1_start, t2))
                {
                    startFound = true;
                    newStart = l1_start;
                }
            }
        }

        if (isPointInTriangle(l1_end, t1) && isPointInTriangle(l1_end, t2))
        {
            endFound = true;
            newEnd = l1_end;
        }

        if (!endFound) {
            if (isPointInTriangle(l1_end, t1))
            {
                if (isPointOnEdge(l1_end, t2))
                {
                    endFound = true;
                    newEnd = l1_end;
                }
            }
        }
    }

    if (secondTime) {
        if (isPointInTriangle(l2_start, t1) && isPointInTriangle(l2_start, t2))
        {
            if (!startFound) {
                startFound = true;
                newStart = l2_start;
            }
            else {
                if (!isTwoPointSame(l1_start, l2_start)) {
                    endFound = true;
                    newEnd = l2_start;
                }
            }
        }

        if (!startFound) {
            if (isPointInTriangle(l2_start, t1))
            {
                if (isPointOnEdge(l2_start, t2))
                {
                    startFound = true;
                    newStart = l2_start;
                }
            }
        }

        if (isPointInTriangle(l2_end, t1) && isPointInTriangle(l2_end, t2))
        {
            if (!endFound) {
                endFound = true;
                newEnd = l2_end;
            }
            else {
                if (!isTwoPointSame(l1_end, l2_end)) {
                    newStart = l2_end;
                    startFound = true;

                }
            
            }
        }

        if (!endFound) {
            if (isPointInTriangle(l2_end, t1))
            {
                if (isPointOnEdge(l2_end, t2))
                {
                    endFound = true;
                    newEnd = l2_end;
                }
            }
        }
    }
    

    if (startFound && endFound)
    {
        pt1 = newStart;
        pt2 = newEnd;
        return true;
    }
    return false;
}


GEOALGOEXPORT int PositionOfThePoint(Vector3D origin, Vector3D normal, Vector3D point) {
    // Calculate D using the plane origin and normal vector
    double D = -(normal.x * origin.x +
        normal.y * origin.y +
        normal.z * origin.z);

    // Calculate the value from the plane equation
    double value = normal.x * point.x +
        normal.y * point.y +
        normal.z * point.z + D;

    // Determine the position relative to the plane
    if (value > 0) {
        return 1;
    }
    else if (value < 0) {
        return -1;
    }
    else {
        return 0;
    }
}
