#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Aabb.h"

#include "Octree.h"

class Triangle
{
    quint64 m_id;
    QVector3D m_normal;
    Aabb m_aabb;
    Octree *m_leaf;
    Octree *m_posInLeaf;

    int m_tagFlag;

    static int tagMask;

public:
    Triangle(quint64 id)
        :m_id(id),
          m_normal(),
          m_aabb(),
          m_leaf(NULL),
          m_posInLeaf(NULL),
          m_tagFlag(1)
    {}

    Triangle(const Triangle &in)
        :m_id(in.m_id),
          m_normal(in.m_normal),
          m_aabb(in.m_aabb),
          m_leaf(in.m_leaf),
          m_posInLeaf(in.m_posInLeaf),
          m_tagFlag(in.m_tagFlag)
    {


    }

    void operator=(const Triangle &in)
    {

        m_id = in.m_id;
        m_normal = in.m_normal;
        m_aabb = in.m_aabb;
        m_leaf = in.m_leaf;
        m_posInLeaf = in.m_posInLeaf;
        m_tagFlag = in.m_tagFlag;
    }
};

#endif // TRIANGLE_H
