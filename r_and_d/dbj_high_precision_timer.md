Source : https://www.codeguru.com/cpp/w-p/system/timers/article.php/c5759/Creating-a-HighPrecision-HighResolution-and-Highly-Reliable-Timer-Utilising-Minimal-CPU-Resources.htm

Posted by: Eugene Manko

# Creating a High-Precision, High-Resolution, and Highly Reliable Timer, Utilising Minimal CPU Resources

> Environment: Windows NT/2K/XP, MSVC 6

While the Microsoft Win32 API provides functions for dealing with waitable timer objects that provide clients with very high resolution (100 nSec.), it gives no guarantee as to the actual precision of those timers. Because those timers rely solely on the APC mechanism to deliver notification to the timer callbacks and those callbacks can only be invoked when threads are put into an APC alertable state, precision is very, very poor. This is because it relies on the thread scheduler, which usually preforms thread context switching every 10msecs delays. Unfortunately, this is not easily overcome. To a lesser extent, the same problem exists when using high-precision multimedia timers.

The proposed solution uses critical sections that don't rely on the thread scheduler to perform time slice allocation and thus are enterable as soon as they become available. This solution utilises the multimedia timer's capabilities to create periodic timers with the highest resolution possible (1msec) and a critical section to provide switching between the threads.

Wait happens because the main program thread gets blocked on the timer critical section while the timer performs counting down cycles. Once this happens, the timer thread leaves the critical section, thus allowing the main thread to continue. Once the main thread enters the imer's critical section, it leaves it immediately, thus allowing the timer thread to block it again on the next pass.

Because this solution does not involve context switching, the delay is minimal. The error level of this timer is below 5%. The presented class can be utilised as shown in the Delay() method below. This example uses high-resolution timers to calculate the actual delay and average delay error level. One can call this method n-number of times with different delay values to verify the error level of this timer.

Enjoy!

```cpp

class PreciseTimer
{
public:
   PreciseTimer() : mRes(0), toLeave(false), stopCounter(-1)
   {
      InitializeCriticalSection(&crit);
      mRes = timeSetEvent(1, 0, &TimerProc, (DWORD)this,
                          TIME_PERIODIC);
   }
   virtual ~PreciseTimer()
   {
      mRes = timeKillEvent(mRes);
      DeleteCriticalSection(&crit);
   }
   
   
   void Wait(int timeout)
   {
      if ( timeout )
      {
         stopCounter = timeout;
         toLeave = true;
         
         EnterCriticalSection(&crit);
         LeaveCriticalSection(&crit);
      }
   }
   
   static void CALLBACK TimerProc(UINT uiID, UINT uiMsg, DWORD dwUser, DWORD dw1, DWORD dw2)
   {
      static volatile bool entered = false;
      
      PreciseTimer* pThis = (PreciseTimer*)dwUser;
      if ( pThis )
      {
         if ( !entered && !pThis->toLeave )   
         {
            entered = true;
            EnterCriticalSection(&pThis->crit);
         }
         else if ( pThis->toLeave && pThis->stopCounter == 0 )
                                              
         {
            pThis->toLeave = false;
            entered = false;
            LeaveCriticalSection(&pThis->crit);
         }
         else if ( pThis->stopCounter > 0 )   
            --pThis->stopCounter;
      }
   }
   
private:
   MMRESULT         mRes;
   CRITICAL_SECTION crit;
   volatile bool    toLeave;
   volatile int     stopCounter;
};


void Delay(unsigned int val)
{
   static LARGE_INTEGER freq = {0};
   static double average = 0;
   static int count = 0;

   ++count;
   LARGE_INTEGER iStart, iStop;
   if ( freq.QuadPart == 0 )
      QueryPerformanceFrequency(&freq), freq.QuadPart /= 1000;
                      

   double sleep = 0;
   QueryPerformanceCounter(&iStart);

   timer.Wait(val); 

   QueryPerformanceCounter(&iStop);
   
   sleep = ((double)iStop.QuadPart - (double)iStart.QuadPart) / (double)freq.QuadPart;
   
   average += (val ? 100.0*(sleep-val)/(double)val : 0);
   
   printf("Waited for %6.3f (%ld). error = %5.2f\n", sleep, val, average/(double)count);
}
```