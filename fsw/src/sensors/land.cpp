#include "../pe_app.h"

// required number of samples for sensor
// to initialize
//
#define 	REQ_LAND_INIT_COUNT (1)
#define 	LAND_TIMEOUT   		(1000000) // 1.0 s

void PE::landInit()
{
	// measure
	math::Vector3F y;
	OS_printf("Land init\n");
	if (landMeasure(y) != CFE_SUCCESS)
	{
		m_LandCount = 0;
	}

	// if finished
	if (m_LandCount > REQ_LAND_INIT_COUNT)
	{
		m_LandTimeout = false;

		(void) CFE_EVS_SendEvent(PE_ESTIMATOR_INF_EID, CFE_EVS_INFORMATION, //TODO update to ld eid
								 "Land detecter init");
	}
}

int PE::landMeasure(math::Vector3F &y)
{
	m_TimeLastLand = m_Timestamp;
	y.Zero();
	m_LandCount += 1;
	return CFE_SUCCESS;
}

void PE::landCorrect()
{
//	// measure land
//	Vector<float, n_y_land> y;
//
//	if (landMeasure(y) != OK) { return; }
//
//	// measurement matrix
//	Matrix<float, n_y_land, n_x> C;
//	C.setZero();
//	// y = -(z - tz)
//	C(Y_land_vx, X_vx) = 1;
//	C(Y_land_vy, X_vy) = 1;
//	C(Y_land_agl, X_z) = -1; // measured altitude, negative down dir.
//	C(Y_land_agl, X_tz) = 1; // measured altitude, negative down dir.
//
//	// use parameter covariance
//	SquareMatrix<float, n_y_land> R;
//	R.setZero();
//	R(Y_land_vx, Y_land_vx) = _land_vxy_stddev.get() * _land_vxy_stddev.get();
//	R(Y_land_vy, Y_land_vy) = _land_vxy_stddev.get() * _land_vxy_stddev.get();
//	R(Y_land_agl, Y_land_agl) = _land_z_stddev.get() * _land_z_stddev.get();
//
//	// residual
//	Matrix<float, n_y_land, n_y_land> S_I = inv<float, n_y_land>((C * _P * C.transpose()) + R);
//	Vector<float, n_y_land> r = y - C * _x;
//	_pub_innov.get().hagl_innov = r(Y_land_agl);
//	_pub_innov.get().hagl_innov_var = R(Y_land_agl, Y_land_agl);
//
//	// fault detection
//	float beta = (r.transpose() * (S_I * r))(0, 0);
//
//	// artifically increase beta threshhold to prevent fault during landing
//	float beta_thresh = 1e2f;
//
//	if (beta / BETA_TABLE[n_y_land] > beta_thresh) {
//		if (!(_sensorFault & SENSOR_LAND)) {
//			_sensorFault |= SENSOR_LAND;
//			mavlink_and_console_log_info(&mavlink_log_pub, "[lpe] land fault,  beta %5.2f", double(beta));
//		}
//
//		// abort correction
//		return;
//
//	} else if (_sensorFault & SENSOR_LAND) {
//		_sensorFault &= ~SENSOR_LAND;
//		mavlink_and_console_log_info(&mavlink_log_pub, "[lpe] land OK");
//	}
//
//	// kalman filter correction always for land detector
//	Matrix<float, n_x, n_y_land> K = _P * C.transpose() * S_I;
//	Vector<float, n_x> dx = K * r;
//	_x += dx;
//	_P -= K * C * _P;
}

void PE::landCheckTimeout()
{
	if (m_Timestamp - m_TimeLastLand > LAND_TIMEOUT)
	{
		if (!m_LandTimeout)
		{
			m_LandTimeout = true;
			m_LandCount = 0;
			(void) CFE_EVS_SendEvent(PE_ESTIMATOR_ERR_EID, CFE_EVS_ERROR,//TODO update to ld eid
									 "Land detecter timeout");
		}
	}
}

