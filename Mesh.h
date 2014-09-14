#ifndef MESH_H
#define MESH_H

#include "Vertex.h"
#include "Triangle.h"
#include "Octree.h"
#include <QVector>
#include <QMatrix4x4>

class Mesh
{
    QVector<Vertex> m_vertices;
    QVector<Triangle> m_triangles;

    QVector<float> m_vertexArray;
    QVector<float> m_normalArray;
    QVector<float> m_colorArray;
    QVector<unsigned int> m_indexArray;

    QVector3D m_center;

    Octree *m_octree;

    QMatrix4x4 m_transformation;

    QVector<Octree *> m_leavesUpdate;

    double m_scale;

    static double globalScale = 500.0;
    static int stateMask = 1;

public:
    Mesh();


    /** Return all the triangles linked to a group of vertices */
     QVector<unsigned int> getTrianglesFromVertices(QVector<unsigned int> &iVerts)
     {
       int triangleTagMask = ++Triangle::tagMask;
       QVector<unsigned int> iTris;
       for (int i = 0; i < iVerts.size(); ++i)
       {
         QVector<unsigned int> ringTris = m_vertices[iVerts[i]].m_adjacentTriangles;
         unsigned int nbTris = ringTris.size();
         for (int j = 0; j < nbTris; ++j)
         {
           unsigned int iTri = ringTris[j];
           if (m_triangles[iTri].m_tagFlag != triangleTagMask)
           {
             iTris.push_back(iTri);
             m_triangles[iTri].m_tagFlag = triangleTagMask;
           }
         }
       }
       return iTris;
     }

     /** Return all the triangles linked to a group of vertices */
     QVector<unsigned int> getVerticesFromTriangles(QVector<unsigned int> &iTris)
     {
       int vertexTagMask = ++Vertex::tagMask;
       QVector<unsigned int> iVerts;
       unsigned int nbTris = iTris.size();

       for (var i = 0; i < nbTris; ++i)
       {
         unsigned int ind = iTris[i] * 3;
         unsigned int iVer1 = iAr[ind];
         unsigned int iVer2 = iAr[ind + 1];
         unsigned int iVer3 = iAr[ind + 2];

         if (m_vertices[iVer1].tagFlag != vertexTagMask)
         {
           iVerts.push_back(iVer1);
           m_vertices[iVer1].tagFlag = vertexTagMask;
         }

         if (vertices[iVer2].tagFlag != vertexTagMask)
         {
           iVerts.push_back(iVer2);
           vertices[iVer2].tagFlag = vertexTagMask;
         }

         if (vertices[iVer3].tagFlag != vertexTagMask)
         {
           iVerts.push_back(iVer3);
           vertices[iVer3].tagFlag = vertexTagMask;
         }
       }
       return iVerts;
     }

     /** Get more triangles (n-ring) */
     void expandsTriangles(QVector<unsigned int> &iTris, unsigned int nRing)
     {
       int triangleTagMask = ++Triangle::tagMask;
       unsigned int nbTris = iTris.size();

       for (int i = 0; i < nbTris; ++i)
       {
         m_triangles[iTris[i]].tagFlag = triangleTagMask;
       }

       int iBegin = 0;
       while (nRing)
       {
         --nRing;
         for (i = iBegin; i < nbTris; ++i)
         {
           unsigned int ind = iTris[i] * 3;
           QVector<unsigned int> iTris1 = m_vertices[iAr[ind]].m_adjacentTriangles;
           QVector<unsigned int> iTris2 = vertices[iAr[ind + 1]].m_adjacentTriangles;
           QVector<unsigned int> iTris3 = vertices[iAr[ind + 2]].m_adjacentTriangles;

           unsigned int nbTris1 = iTris1.size();
           unsigned int nbTris2 = iTris2.size();
           unsigned int nbTris3 = iTris3.size();

           for (int j = 0; j < nbTris1; ++j)
           {
             Triangle t1 = m_triangles[iTris1[j]];
             if (t1.m_tagFlag != triangleTagMask)
             {
               t1.m_tagFlag = triangleTagMask;
               iTris.push_back(iTris1[j]);
             }
           }

           for (int j = 0; j < nbTris2; ++j)
           {
             Triangle t2 = m_triangles[iTris2[j]];
             if (t2.m_tagFlag != triangleTagMask)
             {
               t2.m_tagFlag = triangleTagMask;
               iTris.push_back(iTris2[j]);
             }
           }

           for (int j = 0; j < nbTris3; ++j)
           {
             Triangle t3 = m_triangles[iTris3[j]];
             if (t3.m_tagFlag != triangleTagMask)
             {
               t3.m_tagFlag_ = triangleTagMask;
               iTris.push_back(iTris3[j]);
             }
           }
         }
         iBegin = nbTris;
         nbTris = iTris.size();
       }
     }

     /** Get more vertices (n-ring) */
     void expandsVertices (QVector<unsigned int> &iVerts, unsigned int nRing)
     {
       int vertexTagMask = ++Vertex::tagMask;
       int nbVerts = iVerts.size();

       for (int i = 0; i < nbVerts; ++i)
       {
         m_vertices[iVerts[i]].tagFlag_ = vertexTagMask;
       }

       int iBegin = 0;

       while (nRing)
       {
         --nRing;
         for (int i = iBegin; i < nbVerts; ++i)
         {
           QVector<unsigned int> ring = m_vertices[iVerts[i]].m_ringVertices;
           unsigned int nbRing = ring.size();
           for (int j = 0; j < nbRing; ++j)
           {
             Vertex vRing = vertices[ring[j]];
             if (vRing.m_tagFlag != vertexTagMask)
             {
               vRing.m_tagFlag = vertexTagMask;
               iVerts.push_back(ring[j]);
             }
           }
         }
         iBegin = nbVerts;
         nbVerts = iVerts.size();
       }
     }

     /** Compute the vertices around a vertex */
     computeRingVertices: function (iVert)
     {
       var vertexTagMask = ++Vertex.tagMask_;
       var vertices = this.vertices_;
       var iAr = this.indexArray_;
       var vert = vertices[iVert];
       vert.ringVertices_.length = 0;
       var ring = vert.ringVertices_;
       var iTris = vert.tIndices_;
       var nbTris = iTris.length;
       for (var i = 0; i < nbTris; ++i)
       {
         var ind = iTris[i] * 3;
         var iVer1 = iAr[ind],
           iVer2 = iAr[ind + 1],
           iVer3 = iAr[ind + 2];
         if (iVer1 !== iVert && vertices[iVer1].tagFlag_ !== vertexTagMask)
         {
           ring.push(iVer1);
           vertices[iVer1].tagFlag_ = vertexTagMask;
         }
         if (iVer2 !== iVert && vertices[iVer2].tagFlag_ !== vertexTagMask)
         {
           ring.push(iVer2);
           vertices[iVer2].tagFlag_ = vertexTagMask;
         }
         if (iVer3 !== iVert && vertices[iVer3].tagFlag_ !== vertexTagMask)
         {
           ring.push(iVer3);
           vertices[iVer3].tagFlag_ = vertexTagMask;
         }
       }
     },

     /** Move the mesh center to a certain point */
     moveTo: function (destination)
     {
       mat4.translate(this.matTransform_, mat4.create(), vec3.sub(destination, destination, this.center_));
     },

     /** Render the mesh */
     render: function (camera, picking, lineOrigin, lineNormal)
     {
       this.render_.render(camera, picking, lineOrigin, lineNormal);
     },

     /** Initialize the mesh information : center, octree */
     initMesh: function ()
     {
       var vertices = this.vertices_;
       var triangles = this.triangles_;
       var vAr = this.vertexArray_;
       var nbVertices = vertices.length;
       var nbTriangles = triangles.length;

       //ring vertices and mesh main aabb
       var aabb = new Aabb();
       var i = 0,
         j = 0;
       for (i = 0; i < nbVertices; ++i)
       {
         this.computeRingVertices(i);
         j = i * 3;
         aabb.expandsWithPoint(vAr[j], vAr[j + 1], vAr[j + 2]);
       }
       this.center_ = aabb.computeCenter();

       //scale
       var diag = vec3.dist(aabb.min_, aabb.max_);
       this.scale_ = Mesh.globalScale_ / diag;
       var scale = this.scale_;
       for (i = 0; i < nbVertices; ++i)
       {
         j = i * 3;
         vAr[j] *= scale;
         vAr[j + 1] *= scale;
         vAr[j + 2] *= scale;
       }
       mat4.scale(this.matTransform_, this.matTransform_, [scale, scale, scale]);
       vec3.scale(aabb.min_, aabb.min_, scale);
       vec3.scale(aabb.max_, aabb.max_, scale);
       vec3.scale(this.center_, this.center_, scale);

       //root octree bigger than minimum aabb...
       var vecShift = [0.0, 0.0, 0.0];
       vec3.sub(vecShift, aabb.max_, aabb.min_);
       vec3.scale(vecShift, vecShift, 0.2);
       vec3.sub(aabb.min_, aabb.min_, vecShift);
       vec3.add(aabb.max_, aabb.max_, vecShift);
       aabb.enlargeIfFlat(vec3.length(vecShift)); //for plane mesh...

       //triangles' aabb and normal
       this.updateTrianglesAabbAndNormal();

       //vertices normal
       this.updateVerticesNormal();

       //octree construction
       this.computeOctree(aabb);
     },

     /** compute octree */
     computeOctree: function (aabbSplit)
     {
       var triangles = this.triangles_;
       var nbTriangles = triangles.length;
       var trianglesAll = [];
       for (var i = 0; i < nbTriangles; ++i)
         trianglesAll.push(i);
       this.octree_ = new Octree();
       this.octree_.aabbSplit_.copy(aabbSplit);
       this.octree_.build(this, trianglesAll);
     },

     /** Initialize buffers and shadr */
     initRender: function (shaderType, textures, shaders)
     {
       this.render_.initBuffers();
       this.render_.updateShaders(shaderType, textures, shaders);
       this.render_.updateBuffers();
     },

     /** Update the rendering buffers */
     updateBuffers: function ()
     {
       this.render_.updateBuffers();
     },

     /** Update geometry  */
     updateMesh: function (iTris, iVerts)
     {
       this.updateTrianglesAabbAndNormal(iTris);
       this.updateOctree(iTris);
       this.updateVerticesNormal(iVerts);
     },

     /** Update a group of triangles' normal and aabb */
     updateTrianglesAabbAndNormal: function (iTris)
     {
       var triangles = this.triangles_;
       var vAr = this.vertexArray_;
       var iAr = this.indexArray_;

       var present = iTris !== undefined;
       var nbTris = present ? iTris.length : triangles.length;
       for (var i = 0; i < nbTris; ++i)
       {
         var ind = present ? iTris[i] : i;
         var t = triangles[ind];
         ind *= 3;
         var ind1 = iAr[ind] * 3,
           ind2 = iAr[ind + 1] * 3,
           ind3 = iAr[ind + 2] * 3;
         var v1x = vAr[ind1],
           v1y = vAr[ind1 + 1],
           v1z = vAr[ind1 + 2];
         var v2x = vAr[ind2],
           v2y = vAr[ind2 + 1],
           v2z = vAr[ind2 + 2];
         var v3x = vAr[ind3],
           v3y = vAr[ind3 + 1],
           v3z = vAr[ind3 + 2];
         Geometry.triangleNormal(t.normal_, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z);
         Geometry.computeTriangleAabb(t.aabb_, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z);
       }
     },

     /** Update a group of vertices' normal */
     updateVerticesNormal: function (iVerts)
     {
       var vertices = this.vertices_;
       var triangles = this.triangles_;
       var nAr = this.normalArray_;

       var present = iVerts !== undefined;
       var nbTris = present ? iVerts.length : vertices.length;
       for (var i = 0; i < nbTris; ++i)
       {
         var ind = present ? iVerts[i] : i;
         var vert = vertices[ind];
         var iTris = vert.tIndices_;
         var nbTri = iTris.length;
         var nx = 0.0,
           ny = 0.0,
           nz = 0.0;
         for (var j = 0; j < nbTri; ++j)
         {
           var normTri = triangles[iTris[j]].normal_;
           nx += normTri[0];
           ny += normTri[1];
           nz += normTri[2];
         }
         var len = 1.0 / Math.sqrt(nx * nx + ny * ny + nz * nz);
         ind *= 3;
         nAr[ind] = nx * len;
         nAr[ind + 1] = ny * len;
         nAr[ind + 2] = nz * len;
       }
     },

     /**
      * Update Octree
      * For each triangle we check if its position inside the octree has changed
      * if so... we mark this triangle and we remove it from its former cells
      * We push back the marked triangles into the octree
      */
     updateOctree: function (iTris)
     {
       var nbTris = iTris.length;
       var trisToMove = [];
       var triangles = this.triangles_;
       var i = 0;
       var leaf, trisLeaf;
       for (i = 0; i < nbTris; ++i) //recompute position inside the octree
       {
         var t = triangles[iTris[i]];
         leaf = t.leaf_;
         if (!leaf.aabbSplit_.pointInside(t.aabb_.center_))
         {
           trisToMove.push(iTris[i]);
           trisLeaf = leaf.iTris_;
           if (trisLeaf.length > 0) // remove tris from octree cell
           {
             var iTriLast = trisLeaf[trisLeaf.length - 1];
             var iPos = t.posInLeaf_;
             trisLeaf[iPos] = iTriLast;
             triangles[iTriLast].posInLeaf_ = iPos;
             trisLeaf.pop();
           }
         }
         else if (!t.aabb_.isInside(leaf.aabbLoose_))
           leaf.aabbLoose_.expandsWithAabb(t.aabb_);
       }
       var nbTrisToMove = trisToMove.length;
       for (i = 0; i < nbTrisToMove; ++i) //add triangle to the octree
       {
         var tri = triangles[trisToMove[i]];
         if (this.octree_.aabbLoose_.isOutside(tri.aabb_)) //we reconstruct the whole octree, slow... but rare
         {
           var aabb = new Aabb();
           aabb.copy(this.octree_.aabbSplit_);
           var vecShift = [0.0, 0.0, 0.0];
           vec3.scale(vecShift, vec3.sub(vecShift, aabb.max_, aabb.min_), 0.2);
           vec3.sub(aabb.min_, aabb.min_, vecShift);
           vec3.add(aabb.max_, aabb.max_, vecShift);
           var tris = [];
           var nbTriangles = triangles.length;
           for (i = 0; i < nbTriangles; ++i)
             tris.push(i);
           this.octree_ = new Octree();
           this.octree_.aabbSplit_ = aabb;
           this.octree_.build(this, tris);
           this.leavesUpdate_.length = 0;
           break;
         }
         else
         {
           leaf = tri.leaf_;
           this.octree_.addTriangle(tri);
           if (leaf === tri.leaf_) // failed to insert tri in octree
           {
             trisLeaf = leaf.iTris_;
             tri.posInLeaf_ = trisLeaf.length;
             trisLeaf.push(trisToMove[i]);
           }
         }
       }
     },

     /** End of stroke, update octree (cut empty leaves or go deeper if needed) */
     checkLeavesUpdate: function ()
     {
       Utils.tidy(this.leavesUpdate_);
       var leavesUpdate = this.leavesUpdate_;
       var nbLeaves = leavesUpdate.length;
       var cutLeaves = [];
       var octreeMaxTriangles = Octree.maxTriangles_;
       var octreeMaxDepth = Octree.maxDepth_;
       for (var i = 0; i < nbLeaves; ++i)
       {
         var leaf = leavesUpdate[i];
         if (leaf === null)
           break;
         if (!leaf.iTris_.length)
           leaf.checkEmptiness(cutLeaves);
         else if (leaf.iTris_.length > octreeMaxTriangles && leaf.depth_ < octreeMaxDepth)
           leaf.constructCells(this);
       }
       this.leavesUpdate_.length = 0;
     }
};

#endif // MESH_H
