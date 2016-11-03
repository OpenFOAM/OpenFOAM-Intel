#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <numa.h>
#include <tbb/scalable_allocator.h>
#include <dlfcn.h>

static const size_t alignment = 64;
static uintptr_t start_addr = 0;
static uintptr_t end_addr = 0;

static int is_init = 0;
static int pool_init = 0;
static void* (*real_malloc)(size_t) = 0;
static void* (*real_calloc)(size_t, size_t) = 0;
static void* (*real_realloc)(void*, size_t) = 0;
static void (*real_free)(void*) = 0;

static void real_init(void)
{
    real_malloc = (void* (*)(size_t))dlsym(RTLD_NEXT, "malloc");
    real_calloc = (void* (*)(size_t, size_t))dlsym(RTLD_NEXT, "calloc");
    real_realloc = (void* (*)(void*, size_t))dlsym(RTLD_NEXT, "realloc");
    real_free = (void (*)(void*))dlsym(RTLD_NEXT, "free");

    if (0 == real_malloc || 0 == real_calloc || 0 == real_realloc || 0 == real_free)
    {
        fprintf(stderr, "Error in `dlsym`: %s\n", dlerror());
    }
    
    is_init = 1;
}

static void* raw_alloc(intptr_t /*pool_id*/, size_t& size)
{
    if(start_addr != 0)
    {
        printf("HBMPool: Error - raw_alloc called more than once.\n");
    }
    printf("HBMPool: raw_alloc %ld bytes\n", size);
    void* ptr = numa_alloc_onnode(size, 1);
    if(!ptr)
    {
        printf("raw_alloc: Failed to allocate - NULL pointer returned.\n");
    }
    start_addr = reinterpret_cast<uintptr_t>(ptr);
    end_addr = start_addr + size;
    return ptr;
    
}

static int raw_free(intptr_t /*pool_id*/, void* ptr, size_t size)
{
    printf("HBMPool: raw_free %ld bytes\n", size);
    numa_free(ptr, size);
    return 0;
}

namespace { // private namespace

class HBMPool
{
public:
    HBMPool() : pool(0), hbm_fallback_size(0)
    {
        size_t hbmsize = (size_t)(atoi(getenv("HBM_SIZE")?getenv("HBM_SIZE"):"10")) * 1024L * 1024L;
        threshold = atoi(getenv("HBM_THRESHOLD")?getenv("HBM_THRESHOLD"):"64") * 1024;
        printf("HBMPool: hbmsize=%ld, threshold=%ld\n", hbmsize, threshold);
        rml::MemPoolPolicy pol(raw_alloc, raw_free, hbmsize, true, true);
        rml::pool_create_v1(1, &pol, &pool);
        pool_init = 1;
    }

    ~HBMPool()
    {
        //rml::pool_destroy(pool);
        printf("MBytes requested in HBM but allocated from DDR: %ld\n", hbm_fallback_size/1024/1024);
    }

    void* malloc(size_t size)
    {
        void* ptr = 0;
        if(size > threshold)
        {
            ptr = rml::pool_aligned_malloc(pool, size, alignment);
            if (!ptr)
                hbm_fallback_size += size;
        }
        if(!ptr)
        {
            ptr = real_malloc(size);
        }
        return ptr;
    }
    
    void* realloc(void* ptr, size_t size)
    {
        void* newptr = 0;
        if(0 == size) // free memory
        {
            this->free(ptr);
        }
        else if(0 == ptr) // allocate memory
        {
            newptr = this->malloc(size);
        }
        else if(this->in_pool(ptr)) // in HBM
        {
            if(size > threshold) // is new size greater than threshold then try to realloc in pool
            {
                newptr = rml::pool_aligned_realloc(pool, ptr, size, alignment);
            }
            if(!newptr) // if the size is now smaller or the realloc failed, allocate outside
            {
                newptr = real_malloc(size);
                size_t sz = std::min(size, static_cast<size_t>(end_addr - reinterpret_cast<uintptr_t>(ptr)));
                memcpy(newptr, ptr, sz);
                this->free(ptr);
            }
        }
        else // in DDR
        {
            // no portable way to find out the original size so just realloc on ddt
            // this will mean no reallocations on DDR memory will move to HBM
            newptr = real_realloc(ptr, size);
            if (size > threshold)
            {
                void *newptr2 = rml::pool_aligned_malloc(pool, size, alignment);
                if (newptr2) // if promotion to HBM is successful
                {
                    memcpy(newptr2, newptr, size);
                    real_free(newptr);
                    newptr = newptr2;
                }
            }
        }
        return newptr;
    }
    
    void free(void* ptr)
    {
        if(in_pool(ptr))
        {
            rml::pool_free(pool, ptr);
        }
        else
        {
            real_free(ptr);
        }
    }
    
    bool in_pool(void* ptr)
    {
        const uintptr_t addr = reinterpret_cast<uintptr_t>(ptr);
        return addr >= start_addr && addr < end_addr;
    }
    
private:
    rml::MemoryPool* pool;
    size_t threshold;
    size_t hbm_fallback_size;
};

} // end of private namespace

static HBMPool hbm_pool;

extern "C" {

    void* malloc(size_t size)
    {
        if(0 == is_init)
        {
            real_init();
        }
        if(0 != pool_init)
        {
            return hbm_pool.malloc(size);
        }
        return real_malloc(size);
    }
    void* calloc(size_t nmemb, size_t size)
    {
        if(0 == is_init)
        {
            // Workaround for dlsym() which seems to be happy
            // if calloc() returns NULL
            return NULL;
        }
        if(0 != pool_init)
        {
            void *ptr = hbm_pool.malloc(nmemb*size);
            memset(ptr, 0, nmemb*size);
            return ptr;
        }
        return real_calloc(nmemb, size);
    }
    void* realloc(void* ptr, size_t size)
    {
        if(0 == is_init)
        {
            real_init();
        }
        if(0 != pool_init)
        {
            return hbm_pool.realloc(ptr, size);
        }
        return real_realloc(ptr, size);
    }
    
    void free(void* ptr)
    {
        if(0 == is_init)
        {
            real_init();
        }
        if(0 != pool_init)
        {
            hbm_pool.free(ptr);
        }
        else
        {
            real_free(ptr);
        }
    }

} // extern "C"


