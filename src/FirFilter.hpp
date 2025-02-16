#pragma once
#include <array>
#include <concepts>
#include <cstdint>


template <class T, class C, int size>
    requires(std::integral<T> || std::floating_point<T>)
            && (std::integral<C> || std::floating_point<C> && size > 1)
class FIRFilter
{
private:
    std::array<T, size> m_x;                   // Input buffer
    const std::array<C, size> m_coefficients;  // Coefficients
    T m_output = 0;
    int16_t m_index = 0;

public:
    constexpr FIRFilter(const std::array<C, size>& init)
        : m_coefficients(init)
    {
        m_x.fill(0);
    }

    constexpr T getOutput() const { return m_output; }


    constexpr T filter(T input)
    {
        // Convert input to float for processing
        m_x[m_index] = input;

        C output = 0;
        for (int i = 0; i < size; i++)
        {
            output += m_coefficients[i] * m_x[(m_index + i) % size];
        }
        m_output = output;
        m_index = (m_index + 1) % size;
        return m_output;
    }
};

constexpr bool FIRTest()
{
    std::array<float, 4> coefficients = {0.25, 0.25, 0.25, 0.25};
    FIRFilter<float, float, 4> filter(coefficients);


    filter.filter(1.0);
    filter.filter(0.5);
    filter.filter(2.0);
    filter.filter(1.5);

    if (filter.getOutput() != 1.25)
    {
        return false;
    }
    return true;
}

static_assert(FIRTest(), "Moving average test failed");