#ifndef __EIG3D_H__
#define __EIG3D_H__

#pragma once

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include "Eigen30.h"
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/SVD"
#include "micros.h"

using namespace std;
using namespace Eigen;

void eig3D_eig(const vector<Vector3d>* P,
	const vector<Vector3d>* Q,
	Matrix3d * sigma,
	Matrix3d * rRes,
	Vector3d * tRes,
    double &time);

void eig3D_symbolic(const vector<Vector3d>* P,
	const vector<Vector3d>* Q,
	Matrix3d * sigma,
	Matrix3d * rRes,
	Vector3d * tRes,
    double &time);


#endif
