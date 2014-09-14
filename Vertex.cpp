#include "Vertex.h"

int Vertex::m_tagMask = 1;
int Vertex::m_sculptMask = 1;

Vertex::Vertex(quint64 id):
    m_id(id),
    m_adjacentTriangles(),
    m_ringVertices(),
    m_tagFlag(1),
    m_sculptFlag(1),
    m_stateFlag(1)
{

}
