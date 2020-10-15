/******************************************************************************
**        Title: tRefCPtr.h
**  Description: Smart pointer implementation using reference counting.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#ifndef _TREFCPTR_H_
#define _TREFCPTR_H_

#include <cassert>

namespace gravis
{
  namespace priv
  {


    template <class T>
    class tRefCPtr
    {

      public:
        enum allocType
        {
          ALLOC_OBJECT, ALLOC_ARRAY
        };

        explicit
        tRefCPtr (T* tptr=0, allocType alloct=ALLOC_OBJECT, unsigned int c=1) : refc_ptr(0)
        {
          if (tptr != 0) refc_ptr = new RefCounter(tptr, alloct, c);
        }

        tRefCPtr (const tRefCPtr& rcp)
        {
          _acquireCounter(rcp.refc_ptr);
        }

        ~tRefCPtr ()
        {
          _releaseCounter();
        }

        tRefCPtr& operator= (const tRefCPtr& rcp)
        {
          if (this != &rcp)
          {
            _releaseCounter();
            _acquireCounter(rcp.refc_ptr);
          }
          return *this;
        }

        T& operator*  () const
        {
          return *(refc_ptr->tptr);
        }
        T* operator-> () const
        {
          return refc_ptr->tptr;
        }

        bool isNull() const;

      protected:

        struct RefCounter
        {
            RefCounter (T* ptr=0, allocType _alloct=ALLOC_OBJECT, unsigned int c=1)
              : tptr(ptr), alloct(_alloct), counter(c) {}

          private:
            // Reference Counters should not be copied
            RefCounter(const RefCounter& o) : tptr(o.tptr), alloct(o.alloct), counter(o.counter) {}
            // Reference Counters should not be copied
            RefCounter& operator=(const RefCounter& o)
            {
              tptr=o.tptr;
              alloct=o.alloct;
              counter=o.counter;
            }

          public:
            ~RefCounter ()
            {
              assert(counter == 0);

              if (counter == 0)
              {
                if (alloct == ALLOC_OBJECT)
                  delete tptr;
                else
                  delete[] tptr;
                tptr = 0;
              }
            }

            unsigned int addRef  ()
            {
              return ++counter;
            }
            unsigned int freeRef ()
            {
              return --counter;
            }
            unsigned int getRefCounts () const
            {
              return counter;
            }

            T* tptr;
            allocType    alloct;
            unsigned int counter;
        }* refc_ptr;


        void _acquireCounter (RefCounter* rc)
        {
          refc_ptr = rc;
          if (rc != 0) rc->addRef();
        }

        void _releaseCounter ()
        {
          if (refc_ptr != 0)
          {
            if (refc_ptr->freeRef() == 0)
            {
              delete refc_ptr;
              refc_ptr = 0;
            }
          }
        }
    };

    template <class T>
    inline
    bool tRefCPtr<T>::isNull() const
    {
      return refc_ptr == 0;
    }


  } /* Close namespace "priv" */
} /* Close namespace "gravis" */


#endif /* _TREFCPTR_H_ */
