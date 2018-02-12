#include "../pe_app.h"
#include <math/Matrix1F10.hpp>

/* required number of samples for sensor to initialize */
#define REQ_BARO_INIT_COUNT   (100)
/* 0.1 s */
#define BARO_TIMEOUT          (100000)

void PE::baroInit()
{
    /* measure */
	math::Vector1F y;

	if (baroMeasure(y) != CFE_SUCCESS)
	{
		m_BaroStats.reset();
		return;
	}

	// if finished
	if (m_BaroStats.getCount() > REQ_BARO_INIT_COUNT)
	{
		m_BaroAltOrigin = m_BaroStats.getMean()[0];

		(void) CFE_EVS_SendEvent(PE_SENSOR_INF_EID, CFE_EVS_INFORMATION,
								 "Baro initialized. Mean: (%d) Std dev: (%d) cm",
								 (int)m_BaroStats.getMean()[0],
								 (int)(100 * m_BaroStats.getStdDev()[0]));

		m_BaroTimeout = false;

		if (!m_AltOriginInitialized)
		{
			m_AltOriginInitialized = true;
			m_AltOrigin = m_BaroAltOrigin;
		}
	}
}

int32 PE::baroMeasure(math::Vector1F &y)
{
	//measure
	y.Zero();
	y[0] = m_SensorCombinedMsg.BaroAlt;
	m_BaroStats.update(y);
	m_TimeLastBaro = m_Timestamp;
	return CFE_SUCCESS;
}

void PE::baroCorrect()
{
    /* measure */
    math::Vector1F y;

    if (baroMeasure(y) != CFE_SUCCESS) 
    { 
        return;
    }

    /* subtract baro origin alt */
    y[0] -= m_BaroAltOrigin;



    /* baro measurement matrix */
    math::Matrix1F10 C;
    C.Zero();
    /* measured altitude, negative down dir */
    C[Y_baro_z][X_z] = -1;
//
////	Matrix<float, n_y_baro, n_y_baro> R;
////	R.setZero();
////	R(0, 0) = _baro_stddev.get() * _baro_stddev.get();
//	float residual = 0;
//	residual = ConfigTblPtr->BARO_STDDEV * ConfigTblPtr->BARO_STDDEV;
//
//	float si =
//	// residual
//	Matrix<float, n_y_baro, n_y_baro> S_I =
//		inv<float, n_y_baro>((C * _P * C.transpose()) + R);
//	Vector<float, n_y_baro> r = y - (C * _x);
//
//	//float S_I =
//
//	// fault detection
//	float beta = (r.transpose() * (S_I * r))(0, 0);
//
//	if (beta > BETA_TABLE[n_y_baro]) {
//		if (!(_sensorFault & SENSOR_BARO)) {
//			mavlink_log_critical(&mavlink_log_pub, "[lpe] baro fault, r %5.2f m, beta %5.2f",
//					     double(r(0)), double(beta));
//			_sensorFault |= SENSOR_BARO;
//		}
//
//	} else if (_sensorFault & SENSOR_BARO) {
//		_sensorFault &= ~SENSOR_BARO;
//		mavlink_and_console_log_info(&mavlink_log_pub, "[lpe] baro OK");
//	}
//
//	// kalman filter correction always
//	Matrix<float, n_x, n_y_baro> K = _P * C.transpose() * S_I;
//	Vector<float, n_x> dx = K * r;
//	_x += dx;
//	_P -= K * C * _P;
}

void PE::baroCheckTimeout()
{
	if (m_Timestamp - m_TimeLastBaro > BARO_TIMEOUT)
	{
		if (!m_BaroTimeout)
		{
			m_BaroTimeout = true;
			m_BaroStats.reset();
			(void) CFE_EVS_SendEvent(PE_SENSOR_ERR_EID, CFE_EVS_ERROR,
									 "Baro timeout");
		}
	}
}
