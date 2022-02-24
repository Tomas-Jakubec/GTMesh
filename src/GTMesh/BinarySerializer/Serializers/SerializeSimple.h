#ifndef SERIALIZESIMPLE_H
#define SERIALIZESIMPLE_H
#include <type_traits>
#include <vector>
#include <stdexcept>

namespace Impl {
template<typename T>
constexpr bool IsSimpleSerializable_v = std::is_arithmetic<T>::value || std::is_enum<T>::value;
}

/**
 * @brief The SerializeSimple struct serializes simple types as double or int.
 */
struct SerializeSimple {

    template <typename T, typename ..., std::enable_if_t<Impl::IsSimpleSerializable_v<T>, bool> = true>
    static void serialize(std::vector<unsigned char>& dataContainer, const T& data){
        const unsigned char* pData = reinterpret_cast<const unsigned char*>(&data);
        dataContainer.insert(dataContainer.end(), pData, pData + sizeof (T));
    }

    template <typename T, typename ..., std::enable_if_t<Impl::IsSimpleSerializable_v<T>, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data){
        unsigned char *const pData = reinterpret_cast<unsigned char*>(&data);
        for (size_t i = 0; i < sizeof (T); i++) {
            pData[i] = *dataIterator;
            dataIterator++;
        }
    }

    static void serialize(std::vector<unsigned char>& dataContainer, const bool data){
        dataContainer.push_back(data ? 1 : 0);
    }

    template <typename T, typename ..., std::enable_if_t<Impl::IsSimpleSerializable_v<T>, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, bool& data){
        switch (*dataIterator) {
        case 1: data = true;
            break;
        case 0: data = false;
            break;
        default:
            throw std::runtime_error(std::to_string((int)*dataIterator & 0xff) + ": invalid data input for bool 1/0 expected");
        }
    }

};

#endif // SERIALIZESIMPLE_H
