#include <iostream>
#include "Scheduler.h"

void Scheduler::insertEvent(double lambda, const scheduleEvent &e)
{
    schedule.insert(pair<double, scheduleEvent>(lambda, e));
}

pair<double, scheduleEvent> Scheduler::getNextEvent()
{
    pair<double, scheduleEvent> foo;
    multimap<double, scheduleEvent>::iterator iter;
    
    // get an iterator for the first element
    iter = schedule.begin();
    foo = *iter;
    schedule.erase(iter);
    return(foo);
}



void Scheduler::printSchedule(ostream& outStream)
{
    multimap<double, scheduleEvent>::iterator iter;
    // go through all scheduled events
    for(iter = schedule.begin(); iter!=schedule.end(); ++iter)
    {
        outStream << "Lambda: " << iter->first << endl;
        outStream << "Type: " << iter->second.type << " Group 1: " << iter->second.grp1 << " Group2: " << iter->second.grp2 << endl;
    }
    outStream << endl;
}
