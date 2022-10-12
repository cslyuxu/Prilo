#include "BloomFilter.h"
#include "MurmurHash3.h"

//credited to http://blog.michaelschmatz.com/2016/04/11/how-to-write-a-bloom-filter-cpp/
//https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.h

//https://hur.st/bloomfilter/?n=100000&p=1.0E-6&m=&k=13
//https://drewdevault.com/2016/04/12/How-to-write-a-better-bloom-filter-in-C.html

BloomFilter::BloomFilter(uint64_t size, uint8_t numHashes){
    m_bits.resize(size);
    m_numHashes =  numHashes;
}

void BloomFilter::add(const uint8_t *data, std::size_t len) {
  auto hashValues = hash(data, len);

  for (int n = 0; n < m_numHashes; n++) {
      m_bits[nthHash(n, hashValues[0], hashValues[1], m_bits.size())] = true;
  }
}

bool BloomFilter::possiblyContains(const uint8_t *data, std::size_t len) {
  std::array<uint64_t, 2> hashValues = hash(data, len);

  for (int n = 0; n < m_numHashes; n++) {
      if (!m_bits[nthHash(n, hashValues[0], hashValues[1], m_bits.size())]) {
          return false;
      }
  }

  return true;
}

void BloomFilter::loadBloomFilter(const char *data, uint64_t charsize, uint8_t numHashes){
    m_bits.resize(charsize*8);
    m_numHashes = numHashes;
    char ch;
    for(uint64_t i = 0;i < charsize; i++){
        for(int j = 7;j>=0;j--){
            ch = data[i];
            if((ch&1)==1)
                m_bits[i*8+j] = true;
            ch >> 1;
            
        }
        if(data[i]=='1')
            m_bits[i] = true;
    }
}

void BloomFilter::unloadBloomFilter(char *data, uint64_t size){
    uint64_t charsize = size/8;
    if(size%8!=0){
        abort();
        //printf("The BloomFilter size is not the times of 8!\n");
    }
    char ch;
    for(uint64_t i = 0;i < charsize; i++){
        ch = 0;
        for(int j = 0; j< 8;j++){
            ch = ch << 1 | (int) m_bits[i*8+j];
        }
        data[i] = ch;
    }
}

int BloomFilter::getnumHashes(){
    return m_numHashes;
}

uint64_t BloomFilter::getbitsLen(){
    return m_bits.size();
}

std::array<uint64_t, 2> BloomFilter::hash(const uint8_t *data,std::size_t len){
    std::array<uint64_t, 2> hashValue;
    MurmurHash3_x64_128(data, len, 0, hashValue.data());
    return hashValue;

}

inline uint64_t BloomFilter::nthHash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize) {
    return (hashA + n * hashB) % filterSize;
}




