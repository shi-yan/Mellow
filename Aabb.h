#ifndef AABB_H
#define AABB_H

#include <QVector3D>
#include <QtGlobal>
#include <limits>

class Aabb
{
    QVector3D m_min;
    QVector3D m_max;
    QVector3D m_center;

public:
    Aabb()
    {
        m_min = QVector3D(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        m_max = QVector3D(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
        m_center = QVector3D(0.0, 0.0, 0.0);
    }

    Aabb(const Aabb &in)
    {
        m_min = in.m_min;
        m_max = in.m_max;
        m_center = in.m_center;

    }

    void operator=(const Aabb &in)
    {
        m_min = in.m_min;
        m_max = in.m_max;
        m_center = in.m_center;
    }

    void set(QVector3D min, QVector3D max)
    {
        m_min = min;
        m_max = max;
    }

    void set(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
    {
        m_min = QVector3D(xmin, ymin, zmin);
        m_max = QVector3D(xmax, ymax, zmax);
    }

    QVector3D computeCenter()
    {
        return m_center = (m_min + m_max) * 0.5;
    }


    bool isOutside(const Aabb &in)
    {
        if (in.m_min.x() > m_max.x() || in.m_max.x() < m_min.x())
            return true;
        if (in.m_min.y() > m_max.y() || in.m_max.y() < m_min.y())
            return true;
        if (in.m_min.z() > m_max.z() || in.m_max.z() < m_min.z())
            return true;
        return false;
    }



    /** Return true if aabb is inside the box */
    bool isInside(const Aabb &in)
    {

      if (m_min.x() >= in.m_min.x() && m_max.x() <= in.m_max.x() &&
        m_min.y() >= in.m_min.y() && m_max.y() <= in.m_max.y() &&
        m_min.z() >= in.m_min.z() && m_max.z() <= in.m_max.z())
        return true;
      return false;
    }

    /** Return true if vert is inside the aabb */
    bool pointInside(QVector3D vert)
    {

      if (vert.x() <= m_min.x()) return false;
      if (vert.x() > m_max.x()) return false;
      if (vert.y() <= m_min.y()) return false;
      if (vert.y() > m_max.y()) return false;
      if (vert.z() <= m_min.z()) return false;
      if (vert.z() > m_max.z()) return false;
      return true;
    }

    /** Change the size of the aabb to include vert */
    void expandsWithPoint (double vx, double vy, double vz)
    {
      if (vx > m_max[0]) m_max[0] = vx;
      if (vx < m_min[0]) m_min[0] = vx;
      if (vy > m_max[1]) m_max[1] = vy;
      if (vy < m_min[1]) m_min[1] = vy;
      if (vz > m_max[2]) m_max[2] = vz;
      if (vz < m_min[2]) m_min[2] = vz;
    }

    /** Change the size of the aabb to include another aabb */
    void expandsWithAabb (const Aabb &in)
    {


      if (in.m_max.x() > m_max[0]) m_max[0] = in.m_max.x();
      if (in.m_min.x() < m_min[0]) m_min[0] = in.m_min.x();
      if (in.m_max.y() > m_max[1]) m_max[1] = in.m_max.y();
      if (in.m_min.y() < m_min[1]) m_min[1] = in.m_min.y();
      if (in.m_max.z() > m_max[2]) m_max[2] = in.m_max.z();
      if (in.m_min.z() < m_min[2]) m_min[2] = in.m_min.z();
    }

    /** Return true if a ray intersection the box */
    bool intersectRay(QVector3D vert, QVector3D rayInv)
    {

      double t1 = (m_min[0] - vert.x()) * rayInv.x();
        double t2 = (m_max[0] - vert.x()) * rayInv.x();
        double t3 = (m_min[1] - vert.y()) * rayInv.y();
        double t4 = (m_max[1] - vert.y()) * rayInv.y();
        double t5 = (m_min[2] - vert.z()) * rayInv.z();
        double t6 = (m_max[2] - vert.z()) * rayInv.z();

      double tmin = qMax(qMax(qMin(t1, t2), qMin(t3, t4)), qMin(t5, t6));

       double tmax = qMin(qMin(qMax(t1, t2), qMax(t3, t4)), qMax(t5, t6));
      return (tmax >= 0 && tmin < tmax);
    }

    /** Return true if a sphere intersect with the box */
    bool intersectSphere (QVector3D vert, double radiusSquared)
    {

      double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;

      if (m_min[0] > vert.x())
          dx = m_min[0] - vert.x();
      else if (m_max[0] < vert.x())
          dx = m_max[0] - vert.x();
      else
          dx = 0.0;

      if (m_min[1] > vert.y())
          dy = m_min[1] - vert.y();
      else if (m_max[1] < vert.y())
          dy = m_max[1] - vert.y();
      else dy = 0.0;

      if (m_min[2] > vert.z())
          dz = m_min[2] - vert.z();
      else if (m_max[2] < vert.z())
          dz = m_max[2] - vert.z();
      else dz = 0.0;

      return (dx * dx + dy * dy + dz * dz) < radiusSquared;
    }

    /** Check if the aabb is a plane, if so... enlarge it */
    void enlargeIfFlat(double offset)
    {

      if (m_min[0] == m_max[0])
      {
        m_min[0] -= offset;
        m_max[0] += offset;
      }
      if (m_min[1] == m_max[1])
      {
        m_min[1] -= offset;
        m_max[1] += offset;
      }
      if (m_min[2] == m_max[2])
      {
        m_min[2] -= offset;
        m_max[2] += offset;
      }
    }

    /**
     * Because of laziness I approximate the box by a sphere
     * Return 0 if the sphere is below the plane
     * Return 1 if the sphere is above the plane
     * Return 2 if the sphere intersects the plane
     */
    int intersectPlane (QVector3D origin, QVector3D normal)
    {
      QVector3D center = computeCenter();

      double distToPlane = center.distanceToPlane( origin, normal);
      if (distToPlane * distToPlane <  (m_max - m_min).lengthSquared()  * 0.25)
        return 2;
      return distToPlane > 0.0 ? 1 : 0;
    }


};

#endif // AABB_H
