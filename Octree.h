#ifndef OCTREE_H
#define OCTREE_H

#include <QVector>
#include <Aabb.h>

class Octree
{
    Octree *m_parent;
    int m_depth;
    QVector<Octree*> m_children;
    Aabb m_aabbLoose;
    Aabb m_aabbSplit;
    QVector<quint64> m_triangles;

    static unsigned int maxDepth;
    static unsigned int maxTriangles;

public:
    Octree(Octree *parent = NULL, int depth = 0)
        :m_parent(parent),
          m_depth(depth),
          m_children(),
          m_aabbLoose(),
          m_aabbSplit(),
          m_triangles()
    {

    }



};

#endif // OCTREE_H
