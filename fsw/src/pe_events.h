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

#ifndef PE_EVENTS_H
#define PE_EVENTS_H

#ifdef __cplusplus
extern "C" {
#endif

/************************************************************************
** Local Defines
*************************************************************************/


/* Event IDs
 * Conventions: _EID is the event identifier.  _EVT_CNT is the total number of
 * events and should always be last.  Events can be added before _EVT_CNT.
 * For long-term maintenance, consider not removing an event but replacing it
 * with an unused, reserved, enum to preserve the IDs later in the list. */
typedef enum {

/** \brief <tt> Value of zero is reserved, and should not be used. </tt> */
	PE_RESERVED_EID = 0,  /* Do not use this event ID */

/** \brief <tt> 'Initialized. Version \%d.\%d.\%d.\%d' </tt>
**  \event <tt> 'Initialized. Version \%d.\%d.\%d.\%d' </tt>
**  
**  \par Type: INFORMATION
**
**  \par Cause:
**
**  This event message is issued when the PE task has completed
**  initialization.
**
**  The first \c %d field contains the application major version defined
**      in #PE_MAJOR_VERSION.
**  The second \c %d field contains the application minor version defined
**      in #PE_MINOR_VERSION.
**  The third \c %d field contains the application revision number defined
**      in #PE_REVISION.
**  The fourth \c %d field contains the application revision number defined
**      in #PE_MISSION_REV.
*/
	PE_INIT_INF_EID,

/** \brief <tt> 'Recvd NOOP. Version \%d.\%d.\%d.\%d' </tt>
**  \event <tt> 'Recvd NOOP. Version \%d.\%d.\%d.\%d' </tt>
**  
**  \par Type: INFORMATION
**
**  \par Cause:
**
**  This event message is issued when the CFS PE Task receives a NoOp
**  command.
**
**  The first \c %u field contains the application major version defined
**      in #PE_MAJOR_VERSION.
**  The first \c %d field contains the application major version defined
**      in #PE_MAJOR_VERSION.
**  The second \c %d field contains the application minor version defined
**      in #PE_MINOR_VERSION.
**  The third \c %d field contains the application revision number defined
**      in #PE_REVISION.
**  The fourth \c %d field contains the application revision number defined
**      in #PE_MISSION_REV.
*/
	PE_CMD_NOOP_EID,

    /** \brief <tt> '\%s Pipe failed to subscribe to \%s. (0x\%08X)' </tt>
    **  \event <tt> '\%s Pipe failed to subscribe to \%s. (0x\%08X)' </tt>
**  
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the CFS PE Task fails to subscribe
**  to a message at initialization.
**
**  The first \c %s field contains the name of the pipe.
**  The second \c %s field contains the macro name of the message.
**  The 32 bit hexadecimal number is the error code returned by CFE.
**
*/
	PE_SUBSCRIBE_ERR_EID,

/** \brief <tt> 'Failed to create \%s pipe (0x\%08X)' </tt>
**  \event <tt> 'Failed to create \%s pipe (0x\%08X' </tt>
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the CFS PE Task fails to create a
**  CFE Software Bus pipe at initialization.
**
**  The \c %s field contains the name of the pipe.
**  The 32 bit hexadecimal number is the error code returned by CFE.
**
*/
	PE_PIPE_INIT_ERR_EID,

/** \brief <tt> 'Failed to manage Config table (0x\%08X)' </tt>
**  \event <tt> 'Failed to manage Config table (0x\%08X)' </tt>
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  The 32 bit hexadecimal number is the error code returned by CFE
**  #CFE_TBL_Manage function.
**
*/
	PE_CFGTBL_MANAGE_ERR_EID,

/** \brief <tt> 'Failed to get Config table's address (0x\%08X)' </tt>
**  \event <tt> 'Failed to get Config table's address (0x\%08X)' </tt>
**  
**  \par Type: ERROR
**
**  \par Cause:
**
**  The 32 bit hexadecimal number is the error code returned by CFE
**  #CFE_TBL_GetAddress function.
**
*/
	PE_CFGTBL_GETADDR_ERR_EID,

/** \brief <tt> '\%s pipe read error (0x\%08X).' </tt>
**  \event <tt> '\%s pipe read error (0x\%08X).' </tt>
**  
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the CFS PE Task has
**  had an error reading from a pipe.
**
**  The \c %s field contains the name of the pipe.
**  The 32 bit hexadecimal number is the error code returned by CFE.
**
*/
	PE_RCVMSG_ERR_EID,

/** \brief <tt> 'Recvd invalid \%s msgId (0x\%04X)' </tt>
**  \event <tt> 'Recvd invalid \%s msgId (0x\%04X)' </tt>
**  
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the CFS PE Task has received an invalid
**  message ID.
**
**  The \c %s field contains the name of the pipe.
**  The 16 bit hexadecimal number is the actual message ID received.
*/
	PE_MSGID_ERR_EID,

/** \brief <tt> 'Recvd invalid command code (\%u)' </tt>
**  \event <tt> 'Recvd invalid command code (\%u)' </tt>
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the CFS PE Task has received an invalid
**  message ID.
**
**  The %u field contains the actual command code received.
*/
	PE_CC_ERR_EID,

/** \brief <tt> 'Rcvd invalid msgLen: msgId=0x\%08X, cmdCode=\%d, msgLen=\%d, expectedLen=\%d" </tt>
**  \event <tt> 'Rcvd invalid msgLen: msgId=0x\%08X, cmdCode=\%d, msgLen=\%d, expectedLen=\%d" </tt>
**  
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the CFS PETask has
**  received a message with an invalid length.
**
**  The first \c %d field contains the message ID.
**  The second \c %d field contains the command code.
**  The third \c %d field contains the actual length.
**  The fourth \c %d field contains the expected length.
*/
	PE_MSGLEN_ERR_EID,

/** \brief <tt> 'Failed to register config table (0x%08X)' </tt>
**  \event <tt> 'Failed to register config table (0x%08X)' </tt>
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the CFS PE Task fails to
**  register the PE configuration table.
**
**  The 32 bit hexadecimal number is the error code returned by CFE
**  #CFE_TBL_Register function.
**
*/
	PE_CFGTBL_REG_ERR_EID,

/** \brief <tt> 'Failed to load Config Table (0x%08X)' </tt>
**  \event <tt> 'Failed to load Config Table (0x%08X)' </tt>
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the CFS PE Task fails to
**  load the PE configuration table.
**
**  The 32 bit hexadecimal number is the error code returned by CFE
**  #CFE_TBL_Load function.
**
*/
	PE_CFGTBL_LOAD_ERR_EID,


/** \brief Local estimation initialized
**  \event Local estimation initialized
**
**  \par Type: INFORMATION
**
**  \par Cause:
**
**  This event message is issued when the local estimation has received enough data to initialize
**
*/
	PE_LOCAL_ESTIMATOR_INF_EID,

/** \brief Global estimation initialized
**  \event Global estimation initialized
**
**  \par Type: INFORMATION
**
**  \par Cause:
**
**  This event message is issued when the global estimation has received enough data to initialize
**
*/
	PE_GLOBAL_ESTIMATOR_INF_EID,

/** \brief Reinit state covariance. Index (%i, %i) not finite
**  \event Reinit state covariance. Index (%i, %i) not finite
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when one of the estimator validity checks fails.
**  The estimator will reinitialize after this event is broadcast.
**
*/
	PE_ESTIMATOR_ERR_EID,

/** \brief Baro fault, r %5.2f m, beta %5.2f
**  \event Baro fault, r %5.2f m, beta %5.2f
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the baro fails a validity check.
**  Baro will reinitialize after this event is broadcast.
**
*/
	PE_BARO_FAULT_ERR_EID,

/** \brief Baro initialized. Mean: (%d) Std dev: (%d) cm
**  \event Baro initialized. Mean: (%d) Std dev: (%d) cm
**
**  \par Type: INFORMATION
**
**  \par Cause:
**
**  This event message is issued when baro has (re)initialized successfully.
**
*/

    PE_BARO_OK_INF_EID,

/** \brief Baro timeout
**  \event Baro timeout
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the baro data times out and needs to be reset
**
*/
	PE_BARO_TIMEOUT_ERR_EID,

/** \brief GPS timeout
**  \event GPS timeout
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the gps data times out and needs to be reset
**
*/
	PE_GPS_TIMEOUT_ERR_EID,

/** \brief gps fault, %3g %3g %3g %3g %3g %3g
**  \event gps fault, %3g %3g %3g %3g %3g %3g
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the GPS fails a validity check.
**  GPS will reinitialize after this event is broadcast.
**
*/
	PE_GPS_FAULT_ERR_EID,

/** \brief GPS OK
**  \event GPS OK
**
**  \par Type: INFORMATION
**
**  \par Cause:
**
**  This event message is issued when GPS has (re)initialized successfully.
**
*/
    PE_GPS_OK_INF_EID,

/** \brief Land fault, beta %5.2f
**  \event Land fault, beta %5.2f
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the land detector fails a validity check.
**  Land detector will reinitialize after this event is broadcast.
**
*/
	PE_LAND_FAULT_ERR_EID,

/** \brief Land detector initialized
**  \event Land detector initialized
**
**  \par Type: INFORMATION
**
**  \par Cause:
**
**  This event message is issued when land detector has (re)initialized successfully.
**
*/
	PE_LAND_OK_INF_EID,

/** \brief Land detector timeout
**  \event Land detector timeout
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the land detector data times out and needs to be reset
**
*/
	PE_LAND_TIMEOUT_ERR_EID,

/** \brief PE - mutex creation failed (0x%08lX)
**  \event PE - mutex creation failed (0x%08lX)
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when mutex creation failed.
**
*/
	PE_MUTEX_ERR_EID,

/** \brief Local position message data invalid
**  \event Local position message data invalid
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the local position message data fails validity checks.
**
*/
	PE_LOCAL_POS_MSG_ERR_EID,

/** \brief Global position message data invalid
**  \event Global position message data invalid
**
**  \par Type: ERROR
**
**  \par Cause:
**
**  This event message is issued when the global position message data fails validity checks.
**
*/
	PE_GLOBAL_POS_MSG_ERR_EID,

/** \brief <tt> This is a count of all the app events and should not be used. </tt> */
	PE_EVT_CNT
} PE_EventIds_t;


#ifdef __cplusplus
}
#endif

#endif /* PE_EVENTS_H */

/************************/
/*  End of File Comment */
/************************/
    
