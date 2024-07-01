using DataStructures

struct A 
    x::Float64
    y::Float64
end

Base.isless(a::A, b::A) = a.x < b.x


struct focalcompareA <: Base.Order.Ordering end

Base.lt(o::focalcompareA, a::A, b::A) = 
    (a.x, a.y) < (b.x, b.y)


#now make normal heap with A

h = MutableBinaryMinHeap{A}()
h2 = MutableBinaryHeap{A, focalcompareA}() 

a = A(Int(1), Int(2))


##
maximum(results.time_CBS[results.time_CBS .<= 900])
maximum(results.time_ECBS[results.time_ECBS .<=900])