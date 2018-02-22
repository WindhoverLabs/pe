/****************************************************************************
*
*   Copyright (c) 2017 Windhover Labs, L.L.C. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in
*    the documentation and/or other materials provided with the
*    distribution.
* 3. Neither the name Windhover Labs nor the names of its
*    contributors may be used to endorse or promote products derived
*    from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#include "cfe.h"
#include "pe_app.h"
#include "pe_test_utils.h"
#include "pe_app_test.h"
#include "uttest.h"
#include "ut_osapi_stubs.h"
#include "ut_cfe_sb_stubs.h"
#include "ut_cfe_sb_hooks.h"
#include "ut_cfe_es_stubs.h"
#include "ut_cfe_es_hooks.h"
#include "ut_cfe_evs_stubs.h"
#include "ut_cfe_evs_hooks.h"
#include "ut_cfe_time_stubs.h"
#include "ut_cfe_psp_memutils_stubs.h"
#include "ut_cfe_tbl_stubs.h"
#include "ut_cfe_fs_stubs.h"
#include "ut_cfe_time_stubs.h"
#include <float.h>

int32 hookCalledCount = 0;
int16 DESC_SIZE = 80;

/**************************************************************************
 * Tests for PE_InitEvent()
 **************************************************************************/

/**
 * Test PE_InitEvent() with failed CFE_EVS_Register
 */
void Test_PE_InitEvent_Fail_Register(void)
{
    PE oPE;

    /* Set a fail result for EVS */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_EVENTS_SERVICE | CFE_EVS_NOT_IMPLEMENTED;
    int32 expected = CFE_EVS_APP_NOT_REGISTERED;

    Ut_CFE_EVS_SetReturnCode(UT_CFE_EVS_REGISTER_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitEvent();

    /* Verify results */
    UtAssert_True (result == expected, "InitEvent, failed EVS Register");
}


/**************************************************************************
 * Tests for PE_InitPipe()
 **************************************************************************/
/**
 * Test PE_InitPipe(), fail SCH CFE_SB_CreatePipe
 */
void Test_PE_InitPipe_Fail_CreateSCHPipe(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_CREATEPIPE_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail SB create SCH pipe");
}


/**
 * Test PE_InitPipe(), fail CFE_SB_SubscribeEx for wakeup
 */
void Test_PE_InitPipe_Fail_SubscribeWakeup(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_SUBSCRIBEEX_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail CFE_SB_SubscribeEx for wakeup");
}


/**
 * Test PE_InitPipe(), fail CFE_SB_SubscribeEx for sendhk
 */
void Test_PE_InitPipe_Fail_SubscribeSendHK(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_SUBSCRIBEEX_INDEX, expected, 2);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail CFE_SB_SubscribeEx for sendhk");
}


/**
 * Test PE_InitPipe(), fail CMD CFE_SB_CreatePipe
 */
void Test_PE_InitPipe_Fail_CreateCMDPipe(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_CREATEPIPE_INDEX, expected, 2);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail SB create CMD pipe");
}


/**
 * Test PE_InitPipe(), fail CFE_SB_Subscribe for CMD msg
 */
void Test_PE_InitPipe_Fail_SubscribeCMD(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_SUBSCRIBE_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail CFE_SB_Subscribe for CMD");
}


/**
 * Test PE_InitPipe(), fail DATA CFE_SB_CreatePipe
 */
void Test_PE_InitPipe_Fail_CreateDATAPipe(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_CREATEPIPE_INDEX, expected, 2);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail SB create DATA pipe");
}


/**************************************************************************
 * Tests for PE_InitData()
 **************************************************************************/
/**
 * Test PE_InitData()
 */
void Test_PE_InitData(void)
{
    PE oPE;

    /* Set a fail result */
    int32 expected = CFE_SUCCESS;

    /* Execute the function being tested */
    //oPE.InitData(); //TODO: Test

    /* Verify results */
    //UtAssert_True (result == expected, "InitData");
}


/**************************************************************************
 * Tests for PE_InitApp()
 **************************************************************************/
/**
 * Test PE_InitApp(), fail init event
 */
void Test_PE_InitApp_Fail_InitEvent(void)
{
    PE oPE;

    int32 result = CFE_SUCCESS;
    int32 expected = CFE_EVS_APP_NOT_REGISTERED;

    Ut_CFE_EVS_SetReturnCode(UT_CFE_EVS_REGISTER_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, fail init event");
}


/**
 * Test PE_InitApp(), fail init pipe
 */
void Test_PE_InitApp_Fail_InitPipe(void)
{
    PE oPE;

    int32 result = CFE_SUCCESS;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_CREATEPIPE_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, fail init pipe");
}


/**
 * Test PE_InitApp(), fail init data.
 * NOTE: no current way to fail PE_InitData() in default
 */
void Test_PE_InitApp_Fail_InitData(void)
{
    PE oPE;

    int32 result = CFE_SUCCESS;
    int32 expected = CFE_SUCCESS;

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, fail init data");
}


/**
 * Test PE_InitApp(), fail init config table
 */
void Test_PE_InitApp_Fail_InitConfigTbl(void)
{
    PE oPE;

    int32 result = CFE_SUCCESS;
    int32 expected = CFE_TBL_ERR_INVALID_NAME;

    Ut_CFE_TBL_SetReturnCode(UT_CFE_TBL_REGISTER_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, fail init config table");
}


/**
 * Test PE_InitApp(), Nominal
 */
void Test_PE_InitApp_Nominal(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_EXECUTIVE_SERVICE | CFE_ES_ERR_APP_REGISTER;
    int32 expected = CFE_SUCCESS;

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, nominal");
}

/**************************************************************************
 * Tests for PE_AppMain()
 **************************************************************************/
/**
 * Test PE_AppMain(), Fail RegisterApp
 */
void Test_PE_AppMain_Fail_RegisterApp(void)
{
    PE oPE;

    /* fail the register app */
    Ut_CFE_ES_SetReturnCode(UT_CFE_ES_REGISTERAPP_INDEX, CFE_ES_ERR_APP_REGISTER, 1);

    /* Execute the function being tested */
    oPE.AppMain();
}


/**
 * Test PE_AppMain(), Fail InitApp
 */
void Test_PE_AppMain_Fail_InitApp(void)
{
    PE oPE;

    /* fail the register app */
    Ut_CFE_EVS_SetReturnCode(UT_CFE_EVS_REGISTER_INDEX, CFE_EVS_APP_NOT_REGISTERED, 1);

    /* Execute the function being tested */
    oPE.AppMain();
}


/**
 * Test PE_AppMain(), Fail AcquireConfigPtrs
 */
void Test_PE_AppMain_Fail_AcquireConfigPtrs(void)
{
    PE oPE;

    /* fail the register app */
    Ut_CFE_TBL_SetReturnCode(UT_CFE_TBL_GETADDRESS_INDEX, CFE_TBL_ERR_INVALID_HANDLE, 2);

    /* Execute the function being tested */
    oPE.AppMain();
}


/**
 * Test PE_AppMain(), Invalid Schedule Message
 */
void Test_PE_AppMain_InvalidSchMessage(void)
{
    PE oPE;

    /* The following will emulate behavior of receiving a SCH message to send HK */
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_RCVMSG_INDEX, CFE_SUCCESS, 1);
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_GETMSGID_INDEX, 0, 1);

    Ut_CFE_ES_SetReturnCode(UT_CFE_ES_RUNLOOP_INDEX, FALSE, 2);

    /* Execute the function being tested */
    oPE.AppMain();

}


/**
 * Hook to support: PE_AppMain(), Nominal - SendHK
 */
int32 Test_PE_AppMain_Nominal_SendHK_SendMsgHook(CFE_SB_Msg_t *MsgPtr)
{
    /* TODO:  Test the contents of your HK message here. */

    hookCalledCount++;

    return CFE_SUCCESS;
}

/**
 * Test PE_AppMain(), Nominal - SendHK
 */
void Test_PE_AppMain_Nominal_SendHK(void)
{
    PE oPE;

    /* The following will emulate behavior of receiving a SCH message to WAKEUP */
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_RCVMSG_INDEX, CFE_SUCCESS, 1);
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_GETMSGID_INDEX, PE_SEND_HK_MID, 1);

    Ut_CFE_ES_SetReturnCode(UT_CFE_ES_RUNLOOP_INDEX, FALSE, 2);

    /* Used to verify HK was transmitted correctly. */
    hookCalledCount = 0;
    Ut_CFE_ES_SetFunctionHook(UT_CFE_SB_SENDMSG_INDEX, (void*)&Test_PE_AppMain_Nominal_SendHK_SendMsgHook);

    /* Execute the function being tested */
    oPE.AppMain();

    /* Verify results */
    UtAssert_True (hookCalledCount == 1, "AppMain_Nominal_SendHK");

}


/**
 * Test PE_AppMain(), Nominal - Wakeup
 */
void Test_PE_AppMain_Nominal_Wakeup(void)
{
    PE oPE;

    /* The following will emulate behavior of receiving a SCH message to WAKEUP */
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_RCVMSG_INDEX, CFE_SUCCESS, 1);
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_GETMSGID_INDEX, PE_WAKEUP_MID, 1);

    Ut_CFE_ES_SetReturnCode(UT_CFE_ES_RUNLOOP_INDEX, FALSE, 2);

    /* Execute the function being tested */
    //oPE.AppMain();

}

/**
 * Test PE_AppMain(), Nominal - Land
 */
void Test_PE_AppMain_Land_Nominal_1(void)
{
    PE oPE;
	char desc[DESC_SIZE];

    math::Matrix10F10 input = {
		{0.104776, -0.000000, -0.000000, 0.000320, 0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 0.000000},
		{-0.000000, 0.104776, 0.000000, -0.000000, 0.000320, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000},
		{-0.000000, 0.000000, 0.054885, -0.000000, -0.000000, 0.001135, 0.000000, 0.000000, -0.000001, 0.053125},
		{0.000320, -0.000000, -0.000000, 0.002246, -0.000000, -0.000000, -0.000000, 0.000001, 0.000000, 0.000000},
		{0.000000, 0.000320, -0.000000, -0.000000, 0.002246, -0.000000, -0.000001, -0.000000, -0.000000, 0.000000},
		{0.000000, 0.000000, 0.001135, -0.000000, -0.000000, 0.008013, 0.000000, 0.000001, -0.000004, -0.000035},
		{-0.000000, -0.000000, 0.000000, -0.000000, -0.000001, 0.000000, 0.000006, 0.000000, -0.000000, -0.000000},
		{0.000000, -0.000000, 0.000000, 0.000001, -0.000000, 0.000001, 0.000000, 0.000006, -0.000000, -0.000000},
		{0.000000, -0.000000, -0.000001, 0.000000, -0.000000, -0.000004, -0.000000, -0.000000, 0.000006, 0.000000},
		{0.000000, 0.000000, 0.053125, 0.000000, 0.000000, -0.000035, -0.000000, -0.000000, 0.000000, 0.053222}
    };

    math::Vector10F inputVec = {0.002644, 0.002943, 0.078819, 0.011439, 0.012065, 0.013773, -0.000143, -0.000874, 0.002937, 0.072958};

	math::Matrix10F10 output = {
		{0.104755, 0.000000, 0.000000, 0.000169, 0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.104755, 0.000000, -0.000000, 0.000169, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000},
		{0.000000, 0.000000, 0.053761, -0.000000, -0.000000, 0.000388, 0.000000, 0.000000, -0.000000, 0.053186},
		{0.000169, -0.000000, -0.000000, 0.001183, -0.000000, -0.000000, -0.000000, 0.000001, 0.000000, 0.000000},
		{0.000000, 0.000169, -0.000000, -0.000000, 0.001183, -0.000000, -0.000001, -0.000000, -0.000000, 0.000000},
		{0.000000, 0.000000, 0.000388, -0.000000, -0.000000, 0.007517, 0.000000, 0.000001, -0.000004, 0.000006},
		{-0.000000, -0.000000, 0.000000, -0.000000, -0.000001, 0.000000, 0.000006, 0.000000, -0.000000, 0.000000},
		{0.000000, -0.000000, 0.000000, 0.000001, -0.000000, 0.000001, 0.000000, 0.000006, -0.000000, 0.000000},
		{0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000004, -0.000000, -0.000000, 0.000006, -0.000000},
		{0.000000, 0.000000, 0.053186, 0.000000, 0.000000, 0.000006, 0.000000, 0.000000, -0.000000, 0.053218}
	};

	oPE.m_StateCov = input;
	oPE.m_StateVec = inputVec;

    /* Execute the function being tested */
	oPE.landCorrect();

	/* Verify results */
	for(int i = 0; i < 10; i++)
	{
		for(int j = 0; j < 10; j++)
		{
			sprintf(desc, "Nominal land 1: State cov mismatch at index (%i, %i)", i, j);
			UtAssert_True(fabs(oPE.m_StateCov[i][j] - output[i][j]) < 0.1f, desc);
		}
	}
}

/**
 * Test PE_AppMain(), Nominal - Land
 */
void Test_PE_AppMain_Land_Nominal_2(void)
{
    PE oPE;
    char desc[DESC_SIZE];

    math::Matrix10F10 input = {
		{8.004116, 0.000000, 0.000000, 0.000277, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 8.004116, 0.000000, 0.000000, 0.000277, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.000000, 4.002762, 0.000000, 0.000000, 0.004551, 0.000000, 0.000000, 0.000000, 4.000572},
		{0.000277, 0.000000, 0.000000, 0.002195, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.000277, 0.000000, 0.000000, 0.002195, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.000000, 0.004551, 0.000000, 0.000000, 0.012976, 0.000000, 0.000000, 0.000000, 0.002621},
		{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000000, 0.000000},
		{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000000},
		{0.000000, 0.000000, 4.000572, 0.000000, 0.000000, 0.002621, 0.000000, 0.000000, 0.000000, 4.000137}
    };

    math::Vector10F inputVec = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 9.051394, 0.000000, 0.000000, 0.000000, 0.000000};

    math::Matrix10F10 output = {
		{8.004100, 0.000000, 0.000000, 0.000147, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 8.004100, 0.000000, 0.000000, 0.000147, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.000000, 4.000956, 0.000000, 0.000000, 0.002959, 0.000000, 0.000000, 0.000000, 4.000214},
		{0.000147, 0.000000, 0.000000, 0.001169, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.000147, 0.000000, 0.000000, 0.001169, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.000000, 0.002959, 0.000000, 0.000000, 0.011572, 0.000000, 0.000000, 0.000000, 0.002305},
		{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000000, 0.000000},
		{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000000},
		{0.000000, 0.000000, 4.000214, 0.000000, 0.000000, 0.002305, 0.000000, 0.000000, 0.000000, 4.000066}
    };

    oPE.m_StateCov = input;
    oPE.m_StateVec = inputVec;

    /* Execute the function being tested */
	oPE.landCorrect();

	/* Verify results */
	for(int i = 0; i < 10; i++)
	{
		for(int j = 0; j < 10; j++)
		{
			sprintf(desc, "Nominal land 2: State cov mismatch at index (%i, %i)", i, j);
			UtAssert_True(fabs(oPE.m_StateCov[i][j] - output[i][j]) < 0.01f, desc);
		}
	}

}

/**
 * Test PE_AppMain(), Nominal - Land
 */
void Test_PE_AppMain_Land_Nominal_3(void)
{
	PE oPE;
	char desc[DESC_SIZE];

	math::Matrix10F10 input = {
		{0.045317, -0.000000, 0.000000, 0.000351, 0.000000, -0.000000, -0.000000, 0.000001, 0.000000, 0.000000},
		{-0.000000, 0.045317, -0.000000, -0.000000, 0.000351, -0.000000, -0.000001, -0.000000, -0.000000, 0.000000},
		{0.000000, -0.000000, 0.005979, 0.000000, 0.000000, 0.001213, 0.000000, 0.000001, -0.000003, 0.004090},
		{0.000351, -0.000000, 0.000000, 0.002275, -0.000000, 0.000000, -0.000002, 0.000006, 0.000002, -0.000000},
		{0.000000, 0.000351, 0.000000, -0.000000, 0.002275, 0.000000, -0.000006, -0.000002, -0.000001, 0.000000},
		{-0.000000, -0.000000, 0.001213, 0.000000, 0.000000, 0.008182, 0.000001, 0.000006, -0.000022, -0.000002},
		{-0.000000, -0.000001, 0.000000, -0.000002, -0.000006, 0.000001, 0.000030, 0.000000, -0.000000, -0.000000},
		{0.000001, -0.000000, 0.000001, 0.000006, -0.000002, 0.000006, 0.000000, 0.000030, -0.000000, -0.000000},
		{0.000000, -0.000000, -0.000003, 0.000002, -0.000001, -0.000022, -0.000000, -0.000000, 0.000030, 0.000000},
		{0.000000, 0.000000, 0.004090, -0.000000, 0.000000, -0.000002, -0.000000, -0.000000, 0.000000, 0.004098}
	};

	math::Vector10F inputVec = {0.007867, 0.005784, 0.037025, 0.010325, 0.005103, -0.165113, 0.003133, 0.001280, -0.005324, 0.060882};

	math::Matrix10F10 output = {
		{0.045291, 0.000000, -0.000000, 0.000184, 0.000000, -0.000000, -0.000000, 0.000000, 0.000000, 0.000000},
		{0.000000, 0.045291, -0.000000, -0.000000, 0.000184, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000},
		{-0.000000, -0.000000, 0.004703, 0.000000, 0.000000, 0.000392, 0.000000, 0.000000, -0.000001, 0.004095},
		{0.000184, -0.000000, 0.000000, 0.001191, -0.000000, 0.000000, -0.000001, 0.000003, 0.000001, 0.000000},
		{0.000000, 0.000184, 0.000000, -0.000000, 0.001191, 0.000000, -0.000003, -0.000001, -0.000000, 0.000000},
		{-0.000000, -0.000000, 0.000392, 0.000000, 0.000000, 0.007654, 0.000001, 0.000006, -0.000021, 0.000001},
		{-0.000000, -0.000000, 0.000000, -0.000001, -0.000003, 0.000001, 0.000030, 0.000000, -0.000000, 0.000000},
		{0.000000, -0.000000, 0.000000, 0.000003, -0.000001, 0.000006, 0.000000, 0.000030, -0.000000, 0.000000},
		{0.000000, -0.000000, -0.000001, 0.000001, -0.000000, -0.000021, -0.000000, -0.000000, 0.000030, -0.000000},
		{0.000000, -0.000000, 0.004095, 0.000000, 0.000000, 0.000001, 0.000000, 0.000000, -0.000000, 0.004098}
	};

	oPE.m_StateCov = input;
	oPE.m_StateVec = inputVec;

	/* Execute the function being tested */
	oPE.landCorrect();

	/* Verify results */
	for(int i = 0; i < 10; i++)
	{
		for(int j = 0; j < 10; j++)
		{
			sprintf(desc, "Nominal land 2: State cov mismatch at index (%i, %i)", i, j);
			UtAssert_True(fabs(oPE.m_StateCov[i][j] - output[i][j]) < 0.01f, desc);
		}
	}
}


/**
 * Test PE Baro Correct
 */
void Test_PE_Baro_Correct_Nominal(void)
{
    PE oPE;
    /* Setup initial conditions */
    oPE.InitData();
    oPE.m_Params.BARO_STDDEV = 3.0f;

    oPE.m_BaroAltOrigin = 0.2648827434f;
    oPE.m_SensorCombinedMsg.BaroAlt = 0.2870604694f;
    
    oPE.m_StateCov[0][0] = 0.0519067831f;
    oPE.m_StateCov[0][1] = -0.0f;
    oPE.m_StateCov[0][2] = 0.0f;
    oPE.m_StateCov[0][3] = 0.0002360904f;
    oPE.m_StateCov[0][4] = -0.0000000002f;
    oPE.m_StateCov[0][5] = 0.0000000001f;
    oPE.m_StateCov[0][6] = -0.0000000005f;
    oPE.m_StateCov[0][7] = 0.0000001333f;
    oPE.m_StateCov[0][8] = 0.0000000032f;
    oPE.m_StateCov[0][9] = 0.0f;
    
    oPE.m_StateCov[1][0] = -0.0f;
    oPE.m_StateCov[1][1] = 0.0519067831f;
    oPE.m_StateCov[1][2] = 0.0f;
    oPE.m_StateCov[1][3] = 0.0f;
    oPE.m_StateCov[1][4] = 0.0002360904f;
    oPE.m_StateCov[1][5] = 0.0000000002f;
    oPE.m_StateCov[1][6] = -0.0000001333f;
    oPE.m_StateCov[1][7] = -0.0000000006f;
    oPE.m_StateCov[1][8] = 0.0000000039f;
    oPE.m_StateCov[1][9] = 0.0f;
    
    oPE.m_StateCov[2][0] = 0.0f;
    oPE.m_StateCov[2][1] = 0.0f;
    oPE.m_StateCov[2][2] = 0.0354095176f;
    oPE.m_StateCov[2][3] = -0.0f;
    oPE.m_StateCov[2][4] = -0.0f;
    oPE.m_StateCov[2][5] = 0.0007252924f;
    oPE.m_StateCov[2][6] = -0.0000000168f;
    oPE.m_StateCov[2][7] = 0.0000000128f;
    oPE.m_StateCov[2][8] = -0.0000005278;
    oPE.m_StateCov[2][9] = 0.0343233831f;
    
    oPE.m_StateCov[3][0] = 0.0002360904f;
    oPE.m_StateCov[3][1] = 0.0000000002f;
    oPE.m_StateCov[3][2] = -0.0f;
    oPE.m_StateCov[3][3] = 0.0016303370f;
    oPE.m_StateCov[3][4] = -0.0f;
    oPE.m_StateCov[3][5] = -0.0000000001f;
    oPE.m_StateCov[3][6] = -0.0000000145f;
    oPE.m_StateCov[3][7] = 0.0000013262f;
    oPE.m_StateCov[3][8] = 0.0000000328f;
    oPE.m_StateCov[3][9] = 0.0f;

    oPE.m_StateCov[4][0] = -0.0000000002f;
    oPE.m_StateCov[4][1] = 0.0002360904f;
    oPE.m_StateCov[4][2] = -0.0f;
    oPE.m_StateCov[4][3] = -0.0f;
    oPE.m_StateCov[4][4] = 0.0016303370f;
    oPE.m_StateCov[4][5] = -0.0000000004f;
    oPE.m_StateCov[4][6] = -0.0000013259f;
    oPE.m_StateCov[4][7] = -0.0000000155f;
    oPE.m_StateCov[4][8] = 0.0000000424f;
    oPE.m_StateCov[4][9] = 0.0f;

    oPE.m_StateCov[5][0] = 0.0000000001f;
    oPE.m_StateCov[5][1] = 0.0000000002f;
    oPE.m_StateCov[5][2] = 0.0007252924f;
    oPE.m_StateCov[5][3] = -0.0000000001f;
    oPE.m_StateCov[5][4] = -0.0000000004f;
    oPE.m_StateCov[5][5] = 0.0079122623f;
    oPE.m_StateCov[5][6] = -0.0000001861f;
    oPE.m_StateCov[5][7] = 0.0000001420f;
    oPE.m_StateCov[5][8] = -0.0000058463f;
    oPE.m_StateCov[5][9] = 0.0000023954f;

    oPE.m_StateCov[6][0] = -0.0000000005f;
    oPE.m_StateCov[6][1] = -0.0000001333f;
    oPE.m_StateCov[6][2] = -0.0000000168f;
    oPE.m_StateCov[6][3] = -0.0000000145f;
    oPE.m_StateCov[6][4] = -0.0000013259f;
    oPE.m_StateCov[6][5] = -0.0000001861f;
    oPE.m_StateCov[6][6] = 0.0000084227f;
    oPE.m_StateCov[6][7] = -0.0f;
    oPE.m_StateCov[6][8] = 0.0000000001f;
    oPE.m_StateCov[6][9] = -0.0000000001f;

    oPE.m_StateCov[7][0] = 0.0000001333f;
    oPE.m_StateCov[7][1] = -0.0000000006f;
    oPE.m_StateCov[7][2] = 0.0000000128f;
    oPE.m_StateCov[7][3] = 0.0000013262f;
    oPE.m_StateCov[7][4] = -0.0000000155f;
    oPE.m_StateCov[7][5] = 0.0000001420f;
    oPE.m_StateCov[7][6] = -0.0f;
    oPE.m_StateCov[7][7] = 0.0000084227f;
    oPE.m_StateCov[7][8] = -0.0000000001f;
    oPE.m_StateCov[7][9] = 0.0f;

    oPE.m_StateCov[8][0] = 0.0000000032f;
    oPE.m_StateCov[8][1] = 0.0000000039f;
    oPE.m_StateCov[8][2] = -0.0000005278f;
    oPE.m_StateCov[8][3] = 0.0000000328f;
    oPE.m_StateCov[8][4] = 0.0000000424f;
    oPE.m_StateCov[8][5] = -0.0000058463f;
    oPE.m_StateCov[8][6] = 0.0000000001f;
    oPE.m_StateCov[8][7] = -0.0000000001f;
    oPE.m_StateCov[8][8] = 0.0000084271f;
    oPE.m_StateCov[8][9] = -0.0000000007f;

    oPE.m_StateCov[9][0] = 0.0f;
    oPE.m_StateCov[9][1] = 0.0f;
    oPE.m_StateCov[9][2] = 0.0343233831f;
    oPE.m_StateCov[9][3] = 0.0f;
    oPE.m_StateCov[9][4] = 0.0f;
    oPE.m_StateCov[9][5] = 0.0000023954f;
    oPE.m_StateCov[9][6] = -0.0000000001f;
    oPE.m_StateCov[9][7] = 0.0f;
    oPE.m_StateCov[9][8] = -0.0000000007f;
    oPE.m_StateCov[9][9] = 0.0343402736f;

    oPE.m_StateVec[0] = 0.0028564720f;
    oPE.m_StateVec[1] = 0.0077679274f;
    oPE.m_StateVec[2] = -0.0000669884f;
    oPE.m_StateVec[3] = 0.0075708791f;
    oPE.m_StateVec[4] = 0.0291821100f;
    oPE.m_StateVec[5] = -0.0886378586f;
    oPE.m_StateVec[6] = 0.0004460731f;
    oPE.m_StateVec[7] = -0.0001312319f;
    oPE.m_StateVec[8] = -0.0003727875f;
    oPE.m_StateVec[9] = 0.0085765012f;

    oPE.baroCorrect();
    
    /* Execute the function being tested */
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][0], 0.0519067831f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][1], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][2], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][3], 0.0002360904f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][4], -0.0000000002f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][5], 0.0000000001f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][6], -0.0000000005f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][7], 0.0000001333f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][8], 0.0000000032f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[0][9], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][0], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][1], 0.0519067831f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][2], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][3], 0.0000000002f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][4], 0.0002360904f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][5], 0.0000000002f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][6], -0.0000001333f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][7], -0.0000000006f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][8], 0.0000000039f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[1][9], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][0], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][1], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    printf("oPE.m_StateCov[2][2] %0.10f\n", oPE.m_StateCov[2][2]);
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][2], 0.0352707505f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][3], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][4], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][5], 0.0007224500f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][6], -0.0000000167f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][7], 0.0000000128f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][8], -0.0000005258f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    printf("oPE.m_StateCov[2][9] %0.10f\n", oPE.m_StateCov[2][9]);
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[2][9], 0.0341888703f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][0], 0.0002360904f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][1], 0.0000000002f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][2], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][3], 0.0016303370f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][4], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][5], -0.0000000001f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][6], -0.0000000145f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][7], 0.0000013262f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][8], 0.0000000328f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[3][9], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][0], -0.0000000002f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][1], 0.0002360904f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][2], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][3], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][4], 0.0016303370f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][5], -0.0000000004f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][6], -0.0000013259f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][7], -0.0000000155f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][8], 0.0000000424f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[4][9], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][0], 0.0000000001f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][1], 0.0000000002f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][2], 0.0007224500f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][3], -0.0000000001f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][4], -0.0000000004f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][5], 0.0079122037f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][6], -0.0000001861f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][7], 0.0000001420f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][8], -0.0000058462f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[5][9], -0.0000003598f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][0], -0.0000000005f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][1], -0.0000001333f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][2], -0.0000000167f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][3], -0.0000000145f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][4], -0.0000013259f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][5], -0.0000001861f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][6], 0.0000084227f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][7], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][8], 0.0000000001f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[6][9], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][0], 0.0000001333f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][1], -0.0000000006f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][2], 0.0000000128f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][3], 0.0000013262f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][4], -0.0000000155f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][5], 0.0000001420f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][6], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][7], 0.0000084227f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][8], -0.0000000001f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[7][9], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][0], 0.0000000032f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][1], 0.0000000039f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][2], -0.0000005258f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][3], 0.0000000328f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][4], 0.0000000424f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][5], -0.0000058462f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][6], 0.0000000001f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][7], -0.0000000001f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][8], 0.0000084271f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[8][9], 0.0000000013f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][0], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][1], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
            printf("oPE.m_StateCov[9][2] %.10f\n", oPE.m_StateCov[9][2]);
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][2], 0.0341888703f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][3], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][4], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][5], -0.0000003598f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][6], 0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][7], -0.0f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][8], 0.0000000013f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
        printf("oPE.m_StateCov[9][9] %.10f\n", oPE.m_StateCov[9][9]);
    UtAssert_DoubleCmpAbs(oPE.m_StateCov[9][9], 0.0342098884f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

    UtAssert_DoubleCmpAbs(oPE.m_StateVec[0], 0.0028564720f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[1], 0.0077679274f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[2], -0.0001536398f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[3], 0.0075708791f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[4], 0.0291821100f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[5], -0.0886396319f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[6], 0.0004460731f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[7], -0.0001312319f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[8], -0.0003727862f, 0.00001f, "Test_PE_Baro_Correct_Nominal");
    UtAssert_DoubleCmpAbs(oPE.m_StateVec[9], 0.0084925080f, 0.00001f, "Test_PE_Baro_Correct_Nominal");

}


/**
 * Test PE Baro Multiple Matrix Operations
 */
void Test_PE_Baro_Multiple_Matrix_Ops(void)
{
    math::Matrix10F10 matrix10f10(
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f},
            {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f}
    );
    math::Matrix10F1 matrix10f1(
                {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f}
    );
    math::Matrix1F10 matrix1f10(
                {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f}
    );

    /* 1x10 * 10x10 = 
     * 55 110 165 220 275 330 385 440 495 550
     **/
     
    /* 10x1 * 1x10 = 
        	1	2	3	4	5	6	7	8	9	10
        	2	4	6	8	10	12	14	16	18	20
        	3	6	9	12	15	18	21	24	27	30
        	4	8	12	16	20	24	28	32	36	40
        	5	10	15	20	25	30	35	40	45	50
        	6	12	18	24	30	36	42	48	54	60
        	7	14	21	28	35	42	49	56	63	70
        	8	16	24	32	40	48	56	64	72	80
        	9	18	27	36	45	54	63	72	81	90
        	10	20	30	40	50	60	70	80	90	100
     **/
    matrix10f10 = matrix10f10 - (matrix10f1 * (matrix1f10 * matrix10f10));
    
    for(int i = 0; i<10; i++)
    {
        for(int j = 0; j< 10; j++)
        {
            printf("matrix10f10[i][j] %f\n", matrix10f10[i][j]);
        }
    }

    UtAssert_True (matrix10f10[0][0] == 0.0f, "matrix10f10[0][0] == 0.0f");
    UtAssert_True (matrix10f10[0][1] == 0.0f, "matrix10f10[0][1] == 0.0f");
    UtAssert_True (matrix10f10[0][2] == 0.0f, "matrix10f10[0][2] == 0.0f");
    UtAssert_True (matrix10f10[0][3] == 0.0f, "matrix10f10[0][3] == 0.0f");
    UtAssert_True (matrix10f10[0][4] == 0.0f, "matrix10f10[0][4] == 0.0f");
    UtAssert_True (matrix10f10[0][5] == 0.0f, "matrix10f10[0][5] == 0.0f");
    UtAssert_True (matrix10f10[0][6] == 0.0f, "matrix10f10[0][6] == 0.0f");
    UtAssert_True (matrix10f10[0][7] == 0.0f, "matrix10f10[0][7] == 0.0f");
    UtAssert_True (matrix10f10[0][8] == 0.0f, "matrix10f10[0][8] == 0.0f");
    UtAssert_True (matrix10f10[0][9] == 0.0f, "matrix10f10[0][9] == 0.0f");

    UtAssert_True (matrix10f10[1][0] == -1.0f, "matrix10f10[1][0] == -1.0f");
    UtAssert_True (matrix10f10[1][1] == -2.0f, "matrix10f10[1][1] == -2.0f");
    UtAssert_True (matrix10f10[1][2] == -3.0f, "matrix10f10[1][2] == -3.0f");
    UtAssert_True (matrix10f10[1][3] == -4.0f, "matrix10f10[1][3] == -4.0f");
    UtAssert_True (matrix10f10[1][4] == -5.0f, "matrix10f10[1][4] == -5.0f");
    UtAssert_True (matrix10f10[1][5] == -6.0f, "matrix10f10[1][5] == -6.0f");
    UtAssert_True (matrix10f10[1][6] == -7.0f, "matrix10f10[1][6] == -7.0f");
    UtAssert_True (matrix10f10[1][7] == -8.0f, "matrix10f10[1][7] == -8.0f");
    UtAssert_True (matrix10f10[1][8] == -9.0f, "matrix10f10[1][8] == -9.0f");
    UtAssert_True (matrix10f10[1][9] == -10.0f, "matrix10f10[1][9] == -10.0f");

    UtAssert_True (matrix10f10[2][0] == -2.0f, "matrix10f10[2][0] == -2.0f");
    UtAssert_True (matrix10f10[2][1] == -4.0f, "matrix10f10[2][1] == -4.0f");
    UtAssert_True (matrix10f10[2][2] == -6.0f, "matrix10f10[2][2] == -6.0f");
    UtAssert_True (matrix10f10[2][3] == -8.0f, "matrix10f10[2][3] == -8.0f");
    UtAssert_True (matrix10f10[2][4] == -10.0f, "matrix10f10[2][4] == -10.0f");
    UtAssert_True (matrix10f10[2][5] == -12.0f, "matrix10f10[2][5] == -12.0f");
    UtAssert_True (matrix10f10[2][6] == -14.0f, "matrix10f10[2][6] == -14.0f");
    UtAssert_True (matrix10f10[2][7] == -16.0f, "matrix10f10[2][7] == -16.0f");
    UtAssert_True (matrix10f10[2][8] == -18.0f, "matrix10f10[2][8] == -18.0f");
    UtAssert_True (matrix10f10[2][9] == -20.0f, "matrix10f10[2][9] == -20.0f");

    UtAssert_True (matrix10f10[3][0] == -3.0f, "matrix10f10[3][0] == -3.0f");
    UtAssert_True (matrix10f10[3][1] == -6.0f, "matrix10f10[3][1] == -6.0f");
    UtAssert_True (matrix10f10[3][2] == -9.0f, "matrix10f10[3][2] == -9.0f");
    UtAssert_True (matrix10f10[3][3] == -12.0f, "matrix10f10[3][3] == -12.0f");
    UtAssert_True (matrix10f10[3][4] == -15.0f, "matrix10f10[3][4] == -15.0f");
    UtAssert_True (matrix10f10[3][5] == -18.0f, "matrix10f10[3][5] == -18.0f");
    UtAssert_True (matrix10f10[3][6] == -21.0f, "matrix10f10[3][6] == -21.0f");
    UtAssert_True (matrix10f10[3][7] == -24.0f, "matrix10f10[3][7] == -24.0f");
    UtAssert_True (matrix10f10[3][8] == -27.0f, "matrix10f10[3][8] == -27.0f");
    UtAssert_True (matrix10f10[3][9] == -30.0f, "matrix10f10[3][9] == -30.0f");

    UtAssert_True (matrix10f10[4][0] == -4.0f, "matrix10f10[4][0] == -4.0f");
    UtAssert_True (matrix10f10[4][1] == -8.0f, "matrix10f10[4][1] == -8.0f");
    UtAssert_True (matrix10f10[4][2] == -12.0f, "matrix10f10[4][2] == -12.0f");
    UtAssert_True (matrix10f10[4][3] == -16.0f, "matrix10f10[4][3] == -16.0f");
    UtAssert_True (matrix10f10[4][4] == -20.0f, "matrix10f10[4][4] == -20.0f");
    UtAssert_True (matrix10f10[4][5] == -24.0f, "matrix10f10[4][5] == -24.0f");
    UtAssert_True (matrix10f10[4][6] == -28.0f, "matrix10f10[4][6] == -28.0f");
    UtAssert_True (matrix10f10[4][7] == -32.0f, "matrix10f10[4][7] == -32.0f");
    UtAssert_True (matrix10f10[4][8] == -36.0f, "matrix10f10[4][8] == -36.0f");
    UtAssert_True (matrix10f10[4][9] == -40.0f, "matrix10f10[4][9] == -40.0f");

    UtAssert_True (matrix10f10[5][0] == -5.0f, "matrix10f10[5][0] == -5.0f");
    UtAssert_True (matrix10f10[5][1] == -10.0f, "matrix10f10[5][1] == -10.0f");
    UtAssert_True (matrix10f10[5][2] == -15.0f, "matrix10f10[5][2] == -15.0f");
    UtAssert_True (matrix10f10[5][3] == -20.0f, "matrix10f10[5][3] == -20.0f");
    UtAssert_True (matrix10f10[5][4] == -25.0f, "matrix10f10[5][4] == -25.0f");
    UtAssert_True (matrix10f10[5][5] == -30.0f, "matrix10f10[5][5] == -30.0f");
    UtAssert_True (matrix10f10[5][6] == -35.0f, "matrix10f10[5][6] == -35.0f");
    UtAssert_True (matrix10f10[5][7] == -40.0f, "matrix10f10[5][7] == -40.0f");
    UtAssert_True (matrix10f10[5][8] == -45.0f, "matrix10f10[5][8] == -45.0f");
    UtAssert_True (matrix10f10[5][9] == -50.0f, "matrix10f10[5][9] == -50.0f");

    UtAssert_True (matrix10f10[6][0] == -6.0f, "matrix10f10[6][0] == -6.0f");
    UtAssert_True (matrix10f10[6][1] == -12.0f, "matrix10f10[6][1] == -12.0f");
    UtAssert_True (matrix10f10[6][2] == -18.0f, "matrix10f10[6][2] == -18.0f");
    UtAssert_True (matrix10f10[6][3] == -24.0f, "matrix10f10[6][3] == -24.0f");
    UtAssert_True (matrix10f10[6][4] == -30.0f, "matrix10f10[6][4] == -30.0f");
    UtAssert_True (matrix10f10[6][5] == -36.0f, "matrix10f10[6][5] == -36.0f");
    UtAssert_True (matrix10f10[6][6] == -42.0f, "matrix10f10[6][6] == -42.0f");
    UtAssert_True (matrix10f10[6][7] == -48.0f, "matrix10f10[6][7] == -48.0f");
    UtAssert_True (matrix10f10[6][8] == -54.0f, "matrix10f10[6][8] == -54.0f");
    UtAssert_True (matrix10f10[6][9] == -60.0f, "matrix10f10[6][9] == -60.0f");

    UtAssert_True (matrix10f10[7][0] == -7.0f, "matrix10f10[7][0] == -7.0f");
    UtAssert_True (matrix10f10[7][1] == -14.0f, "matrix10f10[7][1] == -14.0f");
    UtAssert_True (matrix10f10[7][2] == -21.0f, "matrix10f10[7][2] == -21.0f");
    UtAssert_True (matrix10f10[7][3] == -28.0f, "matrix10f10[7][3] == -28.0f");
    UtAssert_True (matrix10f10[7][4] == -35.0f, "matrix10f10[7][4] == -35.0f");
    UtAssert_True (matrix10f10[7][5] == -42.0f, "matrix10f10[7][5] == -42.0f");
    UtAssert_True (matrix10f10[7][6] == -49.0f, "matrix10f10[7][6] == -49.0f");
    UtAssert_True (matrix10f10[7][7] == -56.0f, "matrix10f10[7][7] == -56.0f");
    UtAssert_True (matrix10f10[7][8] == -63.0f, "matrix10f10[7][8] == -63.0f");
    UtAssert_True (matrix10f10[7][9] == -70.0f, "matrix10f10[7][9] == -70.0f");

    UtAssert_True (matrix10f10[8][0] == -8.0f, "matrix10f10[8][0] == -8.0f");
    UtAssert_True (matrix10f10[8][1] == -16.0f, "matrix10f10[8][1] == -16.0f");
    UtAssert_True (matrix10f10[8][2] == -24.0f, "matrix10f10[8][2] == -24.0f");
    UtAssert_True (matrix10f10[8][3] == -32.0f, "matrix10f10[8][3] == -32.0f");
    UtAssert_True (matrix10f10[8][4] == -40.0f, "matrix10f10[8][4] == -40.0f");
    UtAssert_True (matrix10f10[8][5] == -48.0f, "matrix10f10[8][5] == -48.0f");
    UtAssert_True (matrix10f10[8][6] == -56.0f, "matrix10f10[8][6] == -56.0f");
    UtAssert_True (matrix10f10[8][7] == -64.0f, "matrix10f10[8][7] == -64.0f");
    UtAssert_True (matrix10f10[8][8] == -72.0f, "matrix10f10[8][8] == -72.0f");
    UtAssert_True (matrix10f10[8][9] == -80.0f, "matrix10f10[8][9] == -80.0f");

    UtAssert_True (matrix10f10[9][0] == -9.0f, "matrix10f10[9][0] == -9.0f");
    UtAssert_True (matrix10f10[9][1] == -18.0f, "matrix10f10[9][1] == -18.0f");
    UtAssert_True (matrix10f10[9][2] == -27.0f, "matrix10f10[9][2] == -27.0f");
    UtAssert_True (matrix10f10[9][3] == -36.0f, "matrix10f10[9][3] == -36.0f");
    UtAssert_True (matrix10f10[9][4] == -45.0f, "matrix10f10[9][4] == -45.0f");
    UtAssert_True (matrix10f10[9][5] == -54.0f, "matrix10f10[9][5] == -54.0f");
    UtAssert_True (matrix10f10[9][6] == -63.0f, "matrix10f10[9][6] == -63.0f");
    UtAssert_True (matrix10f10[9][7] == -72.0f, "matrix10f10[9][7] == -72.0f");
    UtAssert_True (matrix10f10[9][8] == -81.0f, "matrix10f10[9][8] == -81.0f");
    UtAssert_True (matrix10f10[9][9] == -90.0f, "matrix10f10[9][9] == -90.0f");

}


/**************************************************************************
 * Rollup Test Cases
 **************************************************************************/
void PE_App_Test_AddTestCases(void)
{
    UtTest_Add(Test_PE_InitEvent_Fail_Register, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitEvent_Fail_Register");
    UtTest_Add(Test_PE_InitPipe_Fail_CreateSCHPipe, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_CreateSCHPipe");
    UtTest_Add(Test_PE_InitPipe_Fail_SubscribeWakeup, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_SubscribeWakeup");
    UtTest_Add(Test_PE_InitPipe_Fail_SubscribeSendHK, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_SubscribeSendHK");
    UtTest_Add(Test_PE_InitPipe_Fail_CreateCMDPipe, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_CreateCMDPipe");
    UtTest_Add(Test_PE_InitPipe_Fail_SubscribeCMD, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_SubscribeCMD");
    UtTest_Add(Test_PE_InitPipe_Fail_CreateDATAPipe, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_CreateDATAPipe");
    UtTest_Add(Test_PE_InitData, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitData");
    UtTest_Add(Test_PE_InitApp_Fail_InitEvent, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Fail_InitEvent");
    UtTest_Add(Test_PE_InitApp_Fail_InitPipe, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Fail_InitPipe");
    UtTest_Add(Test_PE_InitApp_Fail_InitData, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Fail_InitData");
    UtTest_Add(Test_PE_InitApp_Fail_InitConfigTbl, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Fail_InitConfigTbl");
    UtTest_Add(Test_PE_InitApp_Nominal, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Nominal");
    UtTest_Add(Test_PE_AppMain_Fail_RegisterApp, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Fail_RegisterApp");
    UtTest_Add(Test_PE_AppMain_Fail_InitApp, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Fail_InitApp");
    UtTest_Add(Test_PE_AppMain_Fail_AcquireConfigPtrs, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Fail_AcquireConfigPtrs");
    UtTest_Add(Test_PE_AppMain_InvalidSchMessage, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_InvalidSchMessage");
    UtTest_Add(Test_PE_AppMain_Nominal_SendHK, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Nominal_SendHK");
    UtTest_Add(Test_PE_AppMain_Nominal_Wakeup, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Nominal_Wakeup");
    UtTest_Add(Test_PE_AppMain_Land_Nominal_1, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Land_Nominal_1");
    UtTest_Add(Test_PE_AppMain_Land_Nominal_2, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Land_Nominal_2");
    UtTest_Add(Test_PE_AppMain_Land_Nominal_3, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Land_Nominal_3");
    UtTest_Add(Test_PE_Baro_Correct_Nominal, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_Baro_Correct_Nominal");
    UtTest_Add(Test_PE_Baro_Multiple_Matrix_Ops, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_Baro_Multiple_Matrix_Ops");
//    UtTest_Add(, PE_Test_Setup, PE_Test_TearDown,
//               "");
//    UtTest_Add(, PE_Test_Setup, PE_Test_TearDown,
//               "");
//    UtTest_Add(, PE_Test_Setup, PE_Test_TearDown,
//               "");
}


