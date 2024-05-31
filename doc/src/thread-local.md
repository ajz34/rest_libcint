# Thread Local Buffer

In function of parallel integral (`cint_crafter.rs`), I need to pre-allocate some memory, as thread-local buffers or caches.

This task should be a routine task in C/C++:
```C++
#pragma omp parallel
{
    // thread-local cache
    double* cache = malloc(size_of_cache * sizeof(double));
    #pragma omp for
    for (int i = 0; i < size_of_i; ++i) {
        /* perform computation with cache */
    }
}
```

In rust, surely there exists some way to perform thread-local tasks
([`thread_local!`](https://doc.rust-lang.org/std/macro.thread_local.html)
or [`Mutex`](https://doc.rust-lang.org/std/sync/struct.Mutex.html)).

However, [`thread_local!`](https://doc.rust-lang.org/std/macro.thread_local.html) requires static buffer allocation,
which does not accept dynamic size of vectors. [`Mutex`](https://doc.rust-lang.org/std/sync/struct.Mutex.html) seems
not efficient enough for this specific task; (probably) most time it just busy locking and unlocking buffers.

All these above is to say, it seems to be extremely difficult for rust to have some thread-local cache.

**So I decided to use unsafe rust.**

**I devoted myself to the unleashing indescribable horrors that shatter my psyche and set my mind adrift in the unknowably infinite cosmos. Hooray!**

Seriously, I'm not sure whether this is a best practice. I hope there could be some more safer way to handle
thread-local buffers efficiently.

- First, I need a function, which converts immutable reference to mutable reference, which is highly unsafe:
    ```rust
    unsafe fn cast_mut_slice<T> (slc: &[T]) -> &mut [T] {
        let len = slc.len();
        let ptr = slc.as_ptr() as *mut T;
        let mut mut_vector = std::slice::from_raw_parts_mut(ptr, len);
        return mut_vector;
    }
    ```
    This function is deliberately not public.

- Second, before parallel region, we allocate **immutable** vectors of buffers:
    ```rust
    let thread_cache: Vec<Vec<f64>> = vec![
        vec![0.; cache_size]; rayon::current_num_threads()];
    ```
    You may challenge me, that this code is still far from efficient.
    In many cases, cache or buffer does not required to be zero-initialized. In C/C++, just `malloc` or
    aligned `malloc` is required. So a more proper way can be
    ```rust
    let thread_cache: Vec<Vec<f64>> = vec![{
        let mut cache = Vec::with_capacity(cache_size);
        unsafe { cache.set_len(cache_size); }
        cache
    }; rayon::current_num_threads()];
    ```
    But I'd say that, in the specific task of electronic integral computation, it's not necessary to
    implement such complicated zero-initialization-free codes.
    Cache or buffer size is quite short; the program needs to frequently retrive cache or buffers.
    So I personally prefer zero-initialization.

    And most importantly, we need `thread_cache` to be **immutable**. Clearly it's not true, since
    we need cache to be modified inside parallel threads. However, we have already using rayon, and
    it does not accept closure (`Fn` type) to accept mutable variables. So as a compromise, we
    decieve the rust compiler by declaring (de facto mutable) `thread_cache` as immutable.

- Third, inside parallel region, retrive mutable cache by thread index:
    ```rust
    (0..n_task).into_par_iter().for_each(|n| {
        let thread_index = current_thread_index().unwrap_or(0);
        let mut cache = unsafe { cast_mut_slice(&thread_cache[thread_index]) };
        /* perform computation with cache */
    }
    ```
    And we obtained mutable cache inside rayon!
    We are careful that these unsafe codes should not incur racing or false-sharing.

As final note, in actual electronic integral tensor evaluation, thread-local buffer or cache is not the
most dangerous part. In both `integral_s1_inplace` and `integral_s2ij_inplace`, we directly write results
into output (`out_with_offset` for `s1`, `out` for `s2ij`). These codes are much more dangerous, and are
prone to hard-to-debug errors.
