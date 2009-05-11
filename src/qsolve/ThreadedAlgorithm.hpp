#include <pthread.h>

using namespace _4ti2_;

void
ThreadedAlgorithm::threaded_compute()
{
    //int code = pthread_create(&id, NULL, ThreadedAlgorithm::entry_point, this);
    pthread_create(&id, NULL, ThreadedAlgorithm::entry_point, this);
}

void
ThreadedAlgorithm::wait()
{
    pthread_join(id, NULL);
}

void*
ThreadedAlgorithm::entry_point(void* ptr)
{
    ThreadedAlgorithm* obj = static_cast<ThreadedAlgorithm*>(ptr);
    obj->compute();
    return 0;
}
