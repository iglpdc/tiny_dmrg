// -*- C++ -*-
/***************************************************************************
 * blitz/reduce.h        Reduction operators: sum, mean, min, max,
 *                       minIndex, maxIndex, product, count, any, all
 *
 * $Id: reduce.h,v 1.8 2008/05/02 10:19:54 papadop Exp $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************/

#ifndef BZ_REDUCE_H
#define BZ_REDUCE_H

#ifndef BZ_BLITZ_H
 #include <blitz/blitz.h>
#endif

#ifndef BZ_NUMTRAIT_H
 #include <blitz/numtrait.h>
#endif

#ifndef BZ_NUMINQUIRE_H
 #include <blitz/numinquire.h>
#endif

//  The various reduce classes.
//  The prototype of the reset method is mandated by the class _bz_ReduceReset
//  in file array/reduce.h

BZ_NAMESPACE(blitz)

template<typename P_sourcetype, typename P_resulttype = BZ_SUMTYPE(P_sourcetype)>
class ReduceSum {
public:

    typedef P_sourcetype T_sourcetype;
    typedef P_resulttype T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceSum() { }

    bool operator()(const T_sourcetype& x,const int=0) { 
        sum_ += x; 
        return true;
    }

    T_resulttype result(const int) const { return sum_; }

    void reset() { sum_ = zero(T_resulttype()); }
 
    static const char* name() { return "sum"; }
 
protected:

    T_resulttype sum_;
};

template<typename P_sourcetype, typename P_resulttype = BZ_FLOATTYPE(P_sourcetype)>
class ReduceMean {
public:

    typedef P_sourcetype T_sourcetype;
    typedef P_resulttype T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceMean() { }

    bool operator()(const T_sourcetype& x,const int=0) { 
        sum_ += x; 
        return true;
    }

    T_resulttype result(const int count) const { return sum_ / count; }

    void reset() { sum_ = zero(T_resulttype()); }

    static const char* name() { return "mean"; }

protected:

    T_resulttype sum_;
};

template<typename P_sourcetype>
class ReduceMin {
public:

    typedef P_sourcetype T_sourcetype;
    typedef P_sourcetype T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceMin() { }

    bool operator()(const T_sourcetype& x,const int=0) {
        if (x < min_)
            min_ = x;
        return true;
    }

    T_resulttype result(const int) const { return min_; }

    void reset() { min_ = huge(P_sourcetype()); }

    static const char* name() { return "min"; }

protected:

    T_resulttype min_;
};

template<typename P_sourcetype>
class ReduceMax {
public:

    typedef P_sourcetype T_sourcetype;
    typedef P_sourcetype T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceMax() { }

    bool operator()(const T_sourcetype& x,const int=0) {
        if (x > max_)
            max_ = x;
        return true;
    }

    T_resulttype result(const int) const { return max_; }

    void reset() { max_ = neghuge(P_sourcetype()); }

    static const char* name() { return "max"; }

protected:

    T_resulttype max_;
};

template <typename T>
struct MinMaxValue {
    void operator=(const T& val) { min = max = val; }
    T min;
    T max;
};

template<typename P_sourcetype>
class ReduceMinMax {
public:

    typedef P_sourcetype              T_sourcetype;
    typedef MinMaxValue<P_sourcetype> T_resulttype;
    typedef T_resulttype              T_numtype;

    static const bool needIndex = false, needInit = true;

    ReduceMinMax() { }

    bool operator()(T_sourcetype x,const int=0) {
        if (x > minmax_.max)
            minmax_.max = x;
        else if (x < minmax_.min)
            minmax_.min = x;
        return true;
    }

    T_resulttype result(int) { return minmax_; }

    void reset(P_sourcetype initialValue) { minmax_ = initialValue; }

    static const char* name() { return "minmax"; }

protected:

    T_resulttype minmax_;
};

template<typename P_sourcetype>
class ReduceMinIndex {
public:

    typedef P_sourcetype T_sourcetype;
    typedef int          T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = true, needInit = false;

    ReduceMinIndex() { }

    bool operator()(const T_sourcetype& x,const T_resulttype& index) {
        if (x < min_) {
            min_ = x;
            index_ = index;
        }
        return true;
    }

    T_resulttype result(const int) const { return index_; }

    void reset(const T_resulttype& index) { 
        min_ = huge(T_sourcetype());
        index_ = index;        
    }

    static const char* name() { return "minIndex"; }

protected:

    T_sourcetype min_;
    T_resulttype index_;
};

template<typename P_sourcetype, int N>
class ReduceMinIndexVector {
public:

    typedef P_sourcetype T_sourcetype;
    typedef TinyVector<int,N> T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = true, needInit = false;

    ReduceMinIndexVector() { }

    bool operator()(const T_sourcetype& x, const T_resulttype& index) {
        if (x < min_) {
            min_ = x;
            index_ = index;
        }
        return true;
    }

    T_resulttype result(const int) const { return index_; }

    void reset(const T_resulttype& index) {
        min_ = huge(T_sourcetype());
        index_ = index;
    }

    static const char* name() { return "minIndexVector"; }

protected:

    T_sourcetype min_;
    T_resulttype index_;
};

template<typename P_sourcetype>
class ReduceMaxIndex {
public:

    typedef P_sourcetype T_sourcetype;
    typedef int          T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = true, needInit = false;

    ReduceMaxIndex() { }

    bool operator()(const T_sourcetype& x,const T_resulttype& index) {
        if (x > max_) {
            max_ = x;
            index_ = index;
        }
        return true;
    }

    T_resulttype result(int) const { return index_; }

    void reset(const T_resulttype& index) {
        max_ = neghuge(T_sourcetype());
        index_ = index;
    }

    static const char* name() { return "maxIndex"; }

protected:

    T_sourcetype max_;
    T_resulttype index_;
};

template<typename P_sourcetype, int N_rank>
class ReduceMaxIndexVector {
public:

    typedef P_sourcetype T_sourcetype;
    typedef TinyVector<int,N_rank> T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = true, needInit = false;

    ReduceMaxIndexVector() { }

    bool operator()(const T_sourcetype& x, const T_resulttype& index) {
        if (x > max_) {
            max_ = x;
            index_ = index;
        }
        return true;
    }

    T_resulttype result(const int) const { return index_; }

    void reset(const T_resulttype& index) {
        max_ = neghuge(T_sourcetype());
        index_ = index;
    }

    static const char* name() { return "maxIndexVector"; }

protected:

    T_sourcetype max_;
    T_resulttype index_;
};

template<typename P_sourcetype>
class ReduceFirst {
public:

    typedef P_sourcetype T_sourcetype;
    typedef int          T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceFirst() { }

    bool operator()(const T_sourcetype& x,const T_resulttype& index) {
        if (x) {
            index_ = index;
            return false;
        } else
            return true;
    }

    T_resulttype result(const int) const { return index_; }

    void reset() { index_ = tiny(int()); }

    static const char* name() { return "first"; }

protected:

    T_resulttype index_;
};

template<typename P_sourcetype>
class ReduceLast {
public:

    typedef P_sourcetype T_sourcetype;
    typedef int          T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceLast() { }

    bool operator()(const T_sourcetype& x,const T_resulttype& index) {
        if (x) {
            index_ = index;
            return true;
        } else
            return true;
    }

    T_resulttype result(const int) const { return index_; }

    void reset() { index_ = huge(int()); }

    static const char* name() { return "last"; }

protected:

    T_resulttype index_;
};

template<typename P_sourcetype, typename P_resulttype = BZ_SUMTYPE(P_sourcetype)>
class ReduceProduct {
public:

    typedef P_sourcetype T_sourcetype;
    typedef P_resulttype T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceProduct() { }

    bool operator()(const T_sourcetype& x,const int=0) { 
        product_ *= x; 
        return true;
    }

    T_resulttype result(const int) const { return product_; }

    void reset() { product_ = one(T_resulttype()); }

    static const char* name() { return "product"; }

protected:

    T_resulttype product_;
};

template<typename P_sourcetype>
class ReduceCount {
public:

    typedef P_sourcetype T_sourcetype;
    typedef int          T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceCount() { }

    bool operator()(const T_sourcetype& x,const int=0) {
        if (bool(x))
            ++count_;
        return true;
    }

    T_resulttype result(const int) const { return count_; }

    void reset() { count_ = zero(T_resulttype()); }

    static const char* name() { return "count"; }

protected:

    T_resulttype count_;
};

template<typename P_sourcetype>
class ReduceAny {
public:

    typedef P_sourcetype T_sourcetype;
    typedef bool     T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceAny() { }

    bool operator()(const T_sourcetype& x,const int=0) {
        if (bool(x)) {
            any_ = true;
            return false;
        }

        return true;
    }

    T_resulttype result(const int) const { return any_; }

    void reset() { any_ = false; }

    static const char* name() { return "any"; }

protected:

    T_resulttype any_;
};

template<typename P_sourcetype>
class ReduceAll {
public:

    typedef P_sourcetype T_sourcetype;
    typedef bool     T_resulttype;
    typedef T_resulttype T_numtype;

    static const bool needIndex = false, needInit = false;

    ReduceAll() { }

    bool operator()(const T_sourcetype& x,const int=0) {
        if (!bool(x)) {
            all_ = false;
            return false;
        } else
            return true;
    }

    T_resulttype result(const int) const { return all_; }

    void reset() { all_ = true; }

    static const char* name() { return "all"; }

protected:

    T_resulttype all_;
}; 

BZ_NAMESPACE_END

#endif // BZ_REDUCE_H
