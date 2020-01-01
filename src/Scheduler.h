#ifndef _SCHEDULER_
#define _SCHEDULER_
#include <map>

using namespace std;

// define the elements of an event to be saved by the scheduler
struct scheduleEvent {
    char type;
    int grp1;
    int grp2;
};



class Scheduler
{
    multimap<double, scheduleEvent> schedule;

public:
    // insert an event into the scheduler
    void insertEvent(double lambda, const scheduleEvent &e);
    // get the next event out, the return value is the value at which it will occur and the event in a pair
    // also delete the event from the schedule at this point; if schedule empty, this will crash
    pair<double, scheduleEvent> getNextEvent();
    
    // check if the schedule is empty
    inline bool empty(){return(schedule.empty());};
    
    // this function prints the whole schedule on the screen; mainly for diagonostic purposes
    void printSchedule(ostream& outStream);
};

#endif
