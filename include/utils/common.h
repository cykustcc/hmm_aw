//
//  common.h
//  hmm_aw
//
//  Created by Yukun Chen on 4/27/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#ifndef common_h
#define common_h

#ifdef __cplusplus
extern "C" {
#endif
  
#define BILLION  1000000000L
  
// Timing, count in nano seconds.
#if defined (_WIN32)
#include <Windows.h>
  
  static inline double getRealTime() {
    FILETIME tm;
    ULONGLONG t;
#if defined(NTDDI_WIN8) && NTDDI_VERSION >= NTDDI_WIN8
    /* Windows 8, Windows Server 2012 and later. ---------------- */
    GetSystemTimePreciseAsFileTime( &tm );
#else
    /* Windows 2000 and later. ---------------------------------- */
    GetSystemTimeAsFileTime( &tm );
#endif
    t = ((ULONGLONG)tm.dwHighDateTime << 32) | (ULONGLONG)tm.dwLowDateTime;
    return (double) t / (double) BILLION;
  }
  
  
#else //_WIN32
#include <time.h>
  
#ifdef __MACH__
#include <sys/time.h>
#include <AvailabilityMacros.h>
#include <mach/clock.h>
#include <mach/mach.h>
  //clock_gettime is not implemented on OSX prior to 10.12
#ifndef MAC_OS_X_VERSION_10_12 //Mac Version not 10.12
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
  static inline int clock_gettime(int clk_id, struct timespec* ts) {
    clock_serv_t cclock;
    mach_timespec_t mts;
    clk_id = 0; // something stupid to get ride of warnings
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts->tv_sec = mts.tv_sec;
    ts->tv_nsec = mts.tv_nsec;
    return 0;
  }
#endif // MAC_OS_X_VERSION_10_12
  static inline double getRealTime() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ( ts.tv_sec ) + (double) ( ts.tv_nsec ) / (double) BILLION ;
  }
#endif // __MACH__
  
#endif // _WIN32
  
#ifdef __cplusplus
}
#endif
#endif /* common_h */
