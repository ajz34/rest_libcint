# Thread Local Buffer

In function of parallel integral (`cint_crafter.rs`), I need to pre-allocate some memory, as
thread-local buffers or caches.

This task should be a simple routine task in C/C++ with OpenMP:
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
([`thread_local!`](https://doc.rust-lang.org/std/macro.thread_local.html). However,
[`thread_local!`](https://doc.rust-lang.org/std/macro.thread_local.html) requires static buffer allocation,
which does not accept dynamic size of vectors.

So a better choice could be [`Mutex`](https://doc.rust-lang.org/std/sync/struct.Mutex.html) or
[`RwLock`](https://doc.rust-lang.org/std/sync/struct.RwLock.html).

The basic idea of `Mutex` and `RwLock` is that when data is retrived by thread-1, then thread-n
(except thread-1) could not retrive this data, protecting data safety.
`RwLock` is more flexible than `Mutex`, in that `Mutex` blocks both simultaneously reading and
writing between threads, while `RwLock` only blocks writing.

I guess that
- `Mutex` is more proper for thread-local cache or buffer allocation;
- `RwLock` is more proper for barrier operations, something like `OMP CRITICAL`.

Taking `Mutex` as example, to perform thread-local allocation and usage,

- First, before parallel region, we allocate `Vec<Mutex<Vec<f64>>>` vectors of buffers:
    ```rust
    let thread_cache = (0..rayon::current_num_threads()).map(|n| {
        Mutex::new(vec![0.; cache_size])
    }).collect_vec();
    ```
    You may challenge me, that this code is still far from efficient.
    In many cases, cache or buffer does not required to be zero-initialized. In C/C++, just `malloc` or
    aligned `malloc` is required. So a more proper way can be
    ```rust
    let thread_cache = (0..rayon::current_num_threads()).map(|n| {
        let mut cache = Vec::with_capacity(cache_size);
        unsafe { cache.set_len(cache_size); }
        Mutex::new(cache)
    }).collect_vec();
    ```
    But I'd say that, in the specific task of electronic integral computation, it's not necessary to
    implement such complicated zero-initialization-free codes.
    Cache or buffer size is quite short; the program needs to frequently retrive cache or buffers.
    So I personally prefer zero-initialization.

- Second, inside parallel region, retrive mutable cache by thread index:
    ```rust
    (0..n_task).into_par_iter().for_each(|n| {
        let thread_index = current_thread_index().unwrap_or(0);
        let mut cache = thread_cache[thread_index].lock().unwrap();
        /* perform computation with cache */
    }
    ```
    And we obtained mutable cache inside rayon!

As final note, in actual electronic integral tensor evaluation, we still need unsafe code.
- Function `self.integral_block::<T>` itself is unsafe.
- In both `integral_s1_inplace` and `integral_s2ij_inplace`, we directly write results
    into output (`out_with_offset` for `s1`, `out` for `s2ij`). These codes are much more dangerous,
    and are prone to hard-to-debug errors.
