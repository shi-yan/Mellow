#ifndef VERTEX_H
#define VERTEX_H

#include <QVector>

class Vertex
{
    quint64 m_id;
    QVector<quint64> m_adjacentTriangles;
    QVector<quint64> m_ringVertices;

    int m_tagFlag;
    int m_sculptFlag;
    int m_stateFlag;

    static int m_tagMask;
    static int m_sculptMask;

public:
    Vertex(quint64 id);

    Vertex(const Vertex &in)
        :m_id(in.m_id),
          m_adjacentTriangles(in.m_adjacentTriangles),
          m_ringVertices(in.m_ringVertices),
          m_tagFlag(in.m_tagFlag),
          m_sculptFlag(in.m_sculptFlag),
          m_stateFlag(in.m_stateFlag)
    {
    }

    void operator=(const Vertex &in)
    {
        m_id = in.m_id;
        m_adjacentTriangles = in.m_adjacentTriangles;
        m_ringVertices = in.m_ringVertices;
        m_tagFlag = in.m_tagFlag;
        m_sculptFlag = in.m_sculptFlag;
        m_stateFlag = in.m_stateFlag;
    }

    void replaceTriangles(quint64 oldId, quint64 newId)
    {
        for(int i = 0;i<m_adjacentTriangles.size();++i)
        {

            if (oldId == m_adjacentTriangles[i])
            {
                m_adjacentTriangles[i] = newId;
                return;
            }
        }
    }

    void replaceRingVertex(quint64 oldId, quint64 newId)
    {
        for(int i = 0;i<m_ringVertices.size();++i)
        {
            if (oldId == m_ringVertices[i])
            {
                m_ringVertices[i] = newId;
            }
        }
    }


      /** Remove triangle */
     void removeTriangle(quint64 index)
      {
        for (int i = 0; i < m_adjacentTriangles.size(); ++i)
        {
          if (index == m_adjacentTriangles[i])
          {
            m_adjacentTriangles[i] = m_adjacentTriangles[m_adjacentTriangles.size() - 1];
            m_adjacentTriangles.pop_back();
            return;
          }
        }
      }

      /** Remove ring vertex */
      void removeRingVertex (quint64 index)
      {
        for (int i = 0; i < m_ringVertices.size(); ++i)
        {
          if (index == m_ringVertices[i])
          {
            m_ringVertices[i] = m_ringVertices[m_ringVertices.size() - 1];
            m_ringVertices.pop_back();
            return;
          }
        }
      }
};

#endif // VERTEX_H
