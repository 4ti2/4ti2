namespace _4ti2_ {

template <class T>
inline
VectorDBase<T>::VectorDBase()
{
    //TODO: start = end = 0;
}

template <class T>
inline
VectorDBase<T>::VectorDBase(Size s)
{
    allocate_memory(s);
}

template <class T>
inline
VectorDBase<T>::VectorDBase(Size s, T v)
{
    assert(s >= 0);
    allocate_memory(s);
    assignT(v);
}

template <class T>
inline
VectorDBase<T>::VectorDBase(T* _start, T* _end)
    : start(_start), end(_end)
{
}

template <class T>
inline
void
VectorDBase<T>::allocate_memory(Size s)
{
    assert(s >= 0);
    start = new T[s];
    end = start+s;
}

template <class T>
inline
void
VectorDBase<T>::free_memory()
{
    delete [] start;
}

template <class T>
inline
const T&
VectorDBase<T>::operator[](Index index) const
{
    assert(index >= 0 && index < get_size());
    return start[index];
}

template <class T>
inline
T&
VectorDBase<T>::operator[](Index index)
{
    assert(index >= 0 && index < get_size());
    return start[index];
}

template <class T>
inline
void
VectorDBase<T>::set(Index index, T v)
{
    assert(index >= 0 && index < get_size());
    start[index] = v;
}


template <class T>
inline
const Size
VectorDBase<T>::get_size() const
{
    return (end-start);
}

template <class T>
inline
bool 
VectorDBase<T>::same(const VectorDBase<T>& v1) const
{
    return (start == v1.start);
}

template <class T>
inline
bool 
VectorDBase<T>::operator==(const VectorDBase<T>& v1) const
{
    assert(get_size() == v1.get_size());
    T* i1 = v1.start;
    for (T* i = start; i != end; ++i) {
        if (*i != *i1) { return false; }
        ++i1;
    }
    return true;
}

template <class T>
inline
bool 
VectorDBase<T>::operator!=(const VectorDBase<T>& v1) const
{
    return !(*this == v1);
}

// Lexicographic ordering.
template <class T>
inline
bool
VectorDBase<T>::operator<(const VectorDBase<T>& v1) const
{
    assert(get_size() == v1.get_size());
    T* i = start; T* i1 = v1.start;
    while (i != end && *i == *i1) { ++i; ++i1; }
    if (i != end && *i < *i1) { return true; }
    return false;
}

template <class T>
inline
void
VectorDBase<T>::swap(VectorDBase<T>& v1)
{
    std::swap(v1.start, start);
    std::swap(v1.end, end);
}

template <class T>
inline
void
VectorDBase<T>::addeq(const VectorDBase<T>& v1)
{
    assert(v1.get_size() == get_size());
    T* i1 = v1.start;
    for (T* i = start; i != end; ++i) { *i += *i1; ++i1; }
}

template <class T>
inline
void
VectorDBase<T>::addeq(const VectorDBase<T>& v1, T m)
{
    assert(v1.get_size() == get_size());
    T* i1 = v1.start;
    for (T* i = start; i != end; ++i) { *i += (*i1)*m; ++i1; }
}

template <class T>
inline
void
VectorDBase<T>::subeq(const VectorDBase<T>& v1)
{
    T* i1 = v1.start;
    for (T* i = start; i != end; ++i) { *i -= *i1; ++i1; }
}

template <class T>
inline
void
VectorDBase<T>::diveq(T d)
{
    for (T* i = start; i != end; ++i) { *i /= d; }
}

template <class T>
inline
void
VectorDBase<T>::muleq(T m)
{
    for (T* i = start; i != end; ++i) { *i *= m; }
}

template <class T>
inline
T
VectorDBase<T>::dot(const VectorDBase<T>& v1) const
{
    T r;
    dot(v1, r);
    return r;
}

template <class T>
inline
void
VectorDBase<T>::dot(const VectorDBase<T>& v1, T& r) const
{
    assert(v1.get_size() == get_size());
    r = 0;
    T* i2 = start;
    for (T* i1 = v1.start; i1 != v1.end; ++i1) { r += (*i1)*(*i2); ++i2; }
}

template <class T>
inline
void
VectorDBase<T>::assignT(T v)
{
    for (T* i = start; i != end; ++i) { *i = v; }
}

template <class T> template <class IndexSet>
inline
void
VectorDBase<T>::assignT(T v, const IndexSet& is)
{
    for (typename IndexSet::Iter i = is.begin(); i != is.end(); ++i) {
        start[i] = v;
    } 
}

template <class T>
inline
void
VectorDBase<T>::assign(const VectorDBase<T>& v1)
{
    assert(get_size() == v1.get_size());
    T* i1 = v1.start;
    for (T* i = start; i != end; ++i) { *i = *i1; ++i1; }
}

template <class T> template <class IndexSet>
inline
void
VectorDBase<T>::assign(const VectorDBase<T>& v1, const IndexSet& is)
{
    assert(v1.get_size() == is.get_size());
    assert(get_size() == is.count());
    T* i = start;
    for (typename IndexSet::Iter i1 = is.begin(); i1 != is.end(); ++i1) {
        *i = v1[i1]; ++i;
    } 
}

template <class T> template <class IndexSet1, class IndexSet2>
inline
void
VectorDBase<T>::assign(const VectorDBase<T>& v1, const IndexSet1& is1, 
        const IndexSet2& is2)
{
    assert(is1.count() == is2.count());
    typename IndexSet2::Iter i2 = is2.begin();
    for (typename IndexSet1::Iter i1 = is1.begin(); i1 != is1.end(); ++i1) {
        start[i2] = v1[i1];
        ++i2;
    } 
}

// r = m1*v1 + m2*v2
template <class T>
inline
void
VectorDBase<T>::add(const VectorDBase<T>& v1, T m1, const VectorDBase<T>& v2, T m2)
{
    assert(v1.get_size() == v2.get_size());
    assert(v1.get_size() == get_size());
    T* i1 = v1.start; T* i2 = v2.start; 
    for (T* i = start; i != end; ++i) { *i = (*i1)*m1 + (*i2)*m2; ++i1; ++i2; }
}

template <class T>
inline
void
VectorDBase<T>::sub(const VectorDBase<T>& v1, const VectorDBase<T>& v2)
{
    assert(v1.get_size() == v2.get_size());
    assert(v1.get_size() == get_size());
    T* i1 = v1.start; T* i2 = v2.start;
    for (T* i = start; i != end; ++i) { *i = *i1 - *i2; ++i1; ++i2; }
}

template <class T>
inline
void
VectorDBase<T>::add(const VectorDBase<T>& v1, const VectorDBase<T>& v2)
{
    assert(v1.get_size() == v2.get_size() && v1.get_size() == get_size());
    T* i1 = v1.start; T* i2 = v2.start; 
    for (T* i = start; i != end; ++i) { *i = *i1 + *i2; ++i1; ++i2; }
}

template <class T>
inline
void
VectorDBase<T>::mul(const VectorDBase<T>& v1, T m)
{
    assert(v1.get_size() == get_size());
    T* i1 = v1.start;
    for (T* i = start; i != end; ++i) { *i = (*i1)*m; ++i1; }
}

template <class T>
inline
void
VectorDBase<T>::div(const VectorDBase<T>& v1, T d)
{
    assert(v1.get_size() == get_size());
    T* i1 = v1.start;
    for (T* i = start; i != end; ++i) { *i = (*i1)/d; ++i1; }
}

template <class T>
inline
void
VectorDBase<T>::normalise()
{
    T* i = start;
    while(i != end && *i == 0) { ++i; }
    if (i == end) return;
    T gcd = *i;
    if (gcd == 1) return;
    ++i;
    while (i != end && *i == 0) { ++i; }
    while (i != end) {
        euclidean(gcd, *i, gcd);
        if (gcd == 1) return;
        ++i;
        while (i != end && *i == 0) { ++i; }
    }
    if (gcd != 1) { diveq(gcd); }
}


template <class T>
inline
VectorD<T>::VectorD()
        : VectorDBase<T>()
{
}

template <class T>
inline
VectorD<T>::VectorD(Size s)
        : VectorDBase<T>(s)
{
}

template <class T>
inline
VectorD<T>::VectorD(const VectorD<T>& v1)
        : VectorDBase<T>(v1.get_size())
{
    VectorDBase<T>::assign(v1);
    *this = v1;
}

template <class T>
inline
VectorD<T>::VectorD(Size s, T v)
        : VectorDBase<T>(s, v)
{
}

template <class T>
inline
VectorD<T>&
VectorD<T>::operator=(const VectorD<T>& v1)
{
    VectorDBase<T>::assign(v1);
    return *this;
}

template <class T>
inline
VectorD<T>::~VectorD()
{
    VectorDBase<T>::free_memory();
}

template <class T>
inline
VectorDLW<T>::VectorDLW()
        : VectorDBase<T>()
{
}

template <class T>
inline
VectorDLW<T>::VectorDLW(const VectorDBase<T>& v1)
        : VectorDBase<T>(v1.start, v1.end)
{
}

template <class T>
inline
VectorDLW<T>::VectorDLW(T* start, T* end)
        : VectorDBase<T>(start, end)
{
}

template <class T>
inline
VectorDLW<T>&
VectorDLW<T>::operator=(const VectorDBase<T>& v1)
{
    VectorDBase<T>::start = v1.start;
    VectorDBase<T>::end = v1.end;
    return *this;
}

} // namespace _4ti2_
