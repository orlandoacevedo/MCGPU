#include <thrust/copy.h>
#include <thrust/remove.h>
#include <thrust/device_ptr.h>
#include <iostream>
#include <iterator>
#include <string>

// this functor returns true if the argument is negative, and false otherwise
struct is_negative
{
    __host__ __device__
    bool operator()(const int x)
    {
        return x < 0;
    }
};


template <typename Iterator>
void print_range(const std::string& name, Iterator first, Iterator last)
{
    typedef typename std::iterator_traits<Iterator>::value_type T;

    std::cout << name << ": ";
    thrust::copy(first, last, std::ostream_iterator<T>(std::cout, " "));  
    std::cout << "\n";
}

int main(void)
{
    // input size
    const int N = 10;

    int h_A[N] = {-2 ,-10, 0, 9, -2, 3, 5, 0, -1, -20};
    int* d_A;

    cudaMalloc(&d_A, sizeof(int) * N);
    cudaMemcpy(d_A, h_A, N * sizeof(int), cudaMemcpyHostToDevice);

    thrust::device_ptr<int> d_thrust_A(d_A);
    print_range("Original", d_thrust_A, d_thrust_A + N);

    // we can also compact sequences with the remove functions, which do the opposite of copy
    thrust::device_ptr<int> d_thrust_end = thrust::remove_if(d_thrust_A, d_thrust_A + N, is_negative());
    
    size_t length = d_thrust_end - d_thrust_A;
    std::cout << std::endl << "Size: " << length << std::endl;
    
    print_range("values", d_thrust_A, d_thrust_end);

    return 0;
}