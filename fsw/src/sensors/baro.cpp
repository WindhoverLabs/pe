#include "../pe_app.h"

void PE::baroInit()
{
    /* Measure */
	math::Vector1F y;

	if (baroMeasure(y) != CFE_SUCCESS)
	{
		m_BaroStats.reset();
		return;
	}

	/* If finished */
	if (m_BaroStats.getCount() > REQ_BARO_INIT_COUNT)
	{
		m_BaroAltOrigin = m_BaroStats.getMean()[0];

		(void) CFE_EVS_SendEvent(PE_BARO_OK_INF_EID, CFE_EVS_INFORMATION,
								 "Baro initialized. Mean: (%d) Std dev: (%d) cm",
								 (int)m_BaroStats.getMean()[0],
								 (int)(100 * m_BaroStats.getStdDev()[0]));

		m_BaroTimeout = FALSE;

		if (!m_AltOriginInitialized)
		{
			m_AltOriginInitialized = TRUE;
			m_AltOrigin = m_BaroAltOrigin;
		}
	}
}


int32 PE::baroMeasure(math::Vector1F &y)
{
	/* Measure */
	y.Zero();
	y[0] = m_SensorCombinedMsg.BaroAlt;
	m_BaroStats.update(y);
	return CFE_SUCCESS;
}


void PE::baroCorrect()
{
    CFE_ES_PerfLogEntry(PE_SENSOR_BARO_PERF_ID);

    if (baroMeasure(m_Baro.y) != CFE_SUCCESS)
    {
        goto end_of_function;
    }

    /* subtract baro origin alt */
    m_Baro.y[0] -= m_BaroAltOrigin;

    /* measured altitude, negative down dir */
    m_Baro.C[Y_baro_z][X_z] = -1.0f;

    m_Baro.R[0][0] = m_Params.BARO_STDDEV * m_Params.BARO_STDDEV;

    /* residual */
    /* ((1x10 * 10x10) * 10x1) + 1x1) */
    m_Baro.S_I = m_Baro.C * m_StateCov * m_Baro.C.Transpose() + m_Baro.R;
    
    /* Take the inverse of a 1x1 matrix (reciprical of the single entry) */
    m_Baro.S_I[0][0] = 1 / m_Baro.S_I[0][0];

    /* Vector1F -  (1x10 * Vector10F) */
    m_Baro.r = m_Baro.y - (m_Baro.C * m_StateVec);

    /* fault detection 1F * 1x1 * 1F */
    m_Baro.beta = m_Baro.r[0] * m_Baro.S_I[0][0] * m_Baro.r[0];

    if (m_Baro.beta > BETA_TABLE[n_y_baro])
    {
        if (!m_BaroFault)
        {
            if(Initialized())
            {
                (void) CFE_EVS_SendEvent(PE_BARO_FAULT_ERR_EID, CFE_EVS_ERROR,
                        "Baro fault, r %5.2f m, beta %5.2f", m_Baro.r[0], m_Baro.beta);
            }
            m_BaroFault = TRUE;
        }

    }
    else if (m_BaroFault)
    {
    	m_BaroFault = FALSE;
        (void) CFE_EVS_SendEvent(PE_BARO_OK_INF_EID, CFE_EVS_INFORMATION,
                "Baro OK");

        m_BaroInitialized = TRUE;
    }

    /* kalman filter correction */

	/* 10x10 * 10x1 * 1x1 */
	m_Baro.K = m_StateCov * m_Baro.C.Transpose() * m_Baro.S_I;

	/* 10x1 * 1x1 */
	m_Baro.temp = m_Baro.K * m_Baro.r;

	m_Baro.dx = m_Baro.temp.ToVector();

	/* 10F + 10F*/
	m_StateVec = m_StateVec + m_Baro.dx;

	/* 10x10 - 10x1 * 1x10 * 10x10 */
	m_StateCov = m_StateCov - m_Baro.K * m_Baro.C * m_StateCov;
end_of_function:

    CFE_ES_PerfLogExit(PE_SENSOR_BARO_PERF_ID);
}


void PE::baroCheckTimeout()
{
    uint64 Timestamp = 0;

	if (m_Timestamp > m_TimeLastBaro)
	{
        Timestamp = m_Timestamp - m_TimeLastBaro;
    }
	else if (m_Timestamp < m_TimeLastBaro)
	{
        Timestamp = m_TimeLastBaro - m_Timestamp;
    }

	if (Timestamp > BARO_TIMEOUT)
	{
		if (!m_BaroTimeout)
		{
			m_BaroTimeout = TRUE;
			m_BaroStats.reset();
			(void) CFE_EVS_SendEvent(PE_BARO_TIMEOUT_ERR_EID, CFE_EVS_ERROR,
									 "Baro timeout: %llu us", Timestamp);
		}
	}
}
