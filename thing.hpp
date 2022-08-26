//how it works? still don't know
//where did I found that?
//https://stackoverflow.com/questions/40546553/call-function-based-on-template-argument-type
//https://stackoverflow.com/a/40547515
//answered Nov 11 '16 at 11:55
//Yakk - Adam Nevraumont

template<class...Fs>
struct overload_t;
// the case where we have a function object:
template<class F>
struct overload_t<F>:F{
  overload_t(F f):F(std::move(f)){}
  using F::operator();
  // boilerplate to ensure these are enabled if possible:
  overload_t(overload_t&&)=default;
  overload_t(overload_t const&)=default;
  overload_t& operator=(overload_t&&)=default;
  overload_t& operator=(overload_t const&)=default;
};
// we cannot inherit from a function pointer.  So
// store one, and write an `operator()` that forwards to it:
template<class R, class...Args>
struct overload_t<R(*)(Args...)>{
  using F=R(*)(Args...);
  F f;
  overload_t(F fin):f(fin){}
  R operator()(Args...args)const{
    return f(std::forward<Args>(args)...);
  }
  overload_t(overload_t&&)=default;
  overload_t(overload_t const&)=default;
  overload_t& operator=(overload_t&&)=default;
  overload_t& operator=(overload_t const&)=default;
};
// the case where we have more than type to overload.
// recursively inherit from the one-arg and the rest-of-arg
// and using operator() to bring both of their () into equal standing:
template<class F0, class...Fs>
struct overload_t<F0,Fs...>:
  overload_t<F0>,
  overload_t<Fs...>
{
  using overload_t<F0>::operator();
  using overload_t<Fs...>::operator();
  overload_t(F0 f0, Fs...fs):
    overload_t<F0>(std::move(f0)),
    overload_t<Fs...>(std::move(fs)...)
  {}
  overload_t(overload_t&&)=default;
  overload_t(overload_t const&)=default;
  overload_t& operator=(overload_t&&)=default;
  overload_t& operator=(overload_t const&)=default;
};
// a helper function to create an overload set without
// having to specify types.  Will be obsolete in C++17:
template<class...Fs>
overload_t<Fs...> overload(Fs...fs){ return {std::move(fs)...};}

//important tprintf:
//awesome example from https://en.cppreference.com/w/cpp/language/parameter_pack
//can use when need a (compile-time and type-safe) format string
//maybe rename it tsprintf

// #include <iostream>
 
// void tprintf(const char* format) // base function
// {
//     std::cout << format;
// }
 
// template<typename T, typename... Targs>
// void tprintf(const char* format, T value, Targs... Fargs) // recursive variadic function
// {
//     for ( ; *format != '\0'; format++ ) {
//         if ( *format == '%' ) {
//            std::cout << value;
//            tprintf(format+1, Fargs...); // recursive call
//            return;
//         }
//         std::cout << *format;
//     }
// }
 
// int main()
// {
//     tprintf("% world% %\n","Hello",'!', 1234, 15.2);
//     return 0;
// }