#-------------------------------------------------
#
# Project created by QtCreator 2014-08-22T16:09:08
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Mellow
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    Vertex.cpp \
    Aabb.cpp \
    Octree.cpp \
    Mesh.cpp \
    Triangle.cpp

HEADERS  += mainwindow.h \
    Vertex.h \
    Aabb.h \
    Octree.h \
    Mesh.h \
    Triangle.h
