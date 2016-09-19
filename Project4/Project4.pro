TEMPLATE = app
CONFIG += console c++11 -o3
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp

# MPI Settings
QMAKE_CXX = mpic++
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(/usr/local/bin/mpic++ --showme:compile)
QMAKE_LFLAGS += $$system(/usr/local/bin/mpic++ --showme:link)
QMAKE_CXXFLAGS += $$system(/usr/local/bin/mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(/usr/local/bin/mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK

HEADERS += \
    lib.h
