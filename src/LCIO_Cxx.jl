__precompile__(false)
module LCIO_Cxx
using Cxx
using StaticArrays

import Base: getindex, start, done, next, length, convert, +

const lcio_path="/afs/desy.de/user/j/jstrube/.julia/v0.6/LCIO/deps/usr/"
addHeaderDir(lcio_path*"include", kind=C_System)
Libdl.dlopen(lcio_path*"lib/liblcio.so", Libdl.RTLD_GLOBAL)
cxxinclude("IO/LCReader.h")
cxx"""
#include "lcio.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/Track.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"
#include "EVENT/LCRelation.h"
#include "EVENT/LCGenericObject.h"
#include <vector>
#include <iostream>
#include <string>
"""

function __init__()
        global reader = icxx"IOIMPL::LCFactory::getInstance()->createLCReader();"
        atexit() do
                icxx"delete $(reader);"
        end
end

immutable Vec
        x::Cdouble
        y::Cdouble
        z::Cdouble
        t::Cdouble
end

+(a::Vec, b::Vec) = Vec(a.x+b.x, a.y+b.y, a.z+b.z, a.t+b.t)

immutable ThreeVec
        x::Cdouble
        y::Cdouble
        z::Cdouble
end

+(a::ThreeVec, b::ThreeVec) = ThreeVec(a.x+b.x, a.y+b.y, a.z+b.z)

immutable CalHit
        x::Cfloat
        y::Cfloat
        z::Cfloat
        E::Cfloat
end

const VectorStringP = cxxt"const std::vector<std::string>*"

start(it::VectorStringP) = 0
next(it::VectorStringP,i) = (it[i], i+1)
done(it::VectorStringP,i) = i >= length(it)
getindex(it::VectorStringP,i) = String( icxx"$(it)->at($i).c_str();" )
length(it::VectorStringP) = icxx"$(it)->size();"

# The iterators use the events as state
# the lcio reader returns NULL at the end of the file
# Because the LCIO getNext function could re-use the memory for the next event,
# the next function should not hold the current and next event at the same time.
# Returning current and reading nextEvent as the new state causes memory corruption
type EventIterator
        current::cxxt"EVENT::LCEvent*"
end
start(it::EventIterator) = C_NULL
next(it::EventIterator, state) = (it.current, C_NULL)
function done(it::EventIterator,state)
        it.current = icxx"$(reader)->readNextEvent();"
        it.current == C_NULL
end
length(it::EventIterator) = icxx"$(reader)->getNumberOfEvents();"

# open file with reader, returns iterator
function open(f::Function, fn::AbstractString)
        icxx"""$(reader)->open($(fn));"""
        try
            f(EventIterator( icxx"(EVENT::LCEvent*)NULL;" ))
        finally
            icxx"""$(reader)->close();"""
        end
        # returns an iterator, initialized with a nullptr
        # the iterator knows about the global reader object
end

# We would like to have a typed collection, but what getCollection returns is unfortunately untyped
# The type is established by reading its name from the collection and mapping it in the LCIOTypemap
immutable LCCollection{T}
        coll::cxxt"EVENT::LCCollection*"
end

const SimCalorimeterHit = cxxt"EVENT::SimCalorimeterHit*"
const TrackerHit = cxxt"EVENT::TrackerHit*"
const SimTrackerHit = cxxt"EVENT::SimTrackerHit*"
const MCParticle = cxxt"EVENT::MCParticle*"
const Track = cxxt"EVENT::Track*"
const LCGenericObject = cxxt"EVENT::LCGenericObject*"
const ReconstructedParticle = cxxt"EVENT::ReconstructedParticle*"

# map from names stored in collection to actual types
LCIOTypemap = Dict(
        "SimCalorimeterHit" => SimCalorimeterHit,
        "TrackerHit" => TrackerHit,
        "SimTrackerHit" => SimTrackerHit,
        "MCParticle" => MCParticle,
        "Track" => Track,
        "LCGenericObject" => LCGenericObject,
        "ReconstructedParticle" => ReconstructedParticle
)

start(it::LCCollection) = 0
done(it::LCCollection, i) = i >= length(it)
next{T}(it::LCCollection{T}, i) = icxx"static_cast<$(T)>($(it.coll)->getElementAt($(i)));", i+1
length(it::LCCollection) = icxx"$(it.coll)->getNumberOfElements();"

function getCollection(event, collectionName)
    collection = icxx"""$(event)->getCollection($(pointer(collectionName)));"""
    collectionType = icxx"$(collection)->getTypeName();"
    return LCCollection{LCIOTypemap[String(collectionType)]}(collection)
end

getCollectionTypeName(collection::LCCollection) = String( icxx"$(collection.coll)->getTypeName().c_str();" )
getCollectionNames(event) = icxx"$(event)->getCollectionNames();"

getEnergy(particle) = icxx"$(particle)->getEnergy();"
function getMomentum(particle)
        p3 = icxx"$(particle)->getMomentum();"
        SVector{3, Float64}(unsafe_load(p3, 1), unsafe_load(p3, 2), unsafe_load(p3, 3))
end
function getPosition(hit)
        pos = icxx"$(hit)->getPosition();"
        SVector{3, Float64}(unsafe_load(pos, 1), unsafe_load(pos, 2), unsafe_load(pos, 3))
end

function getP4(particle)
        p3 = icxx"$(particle)->getMomentum();"
        e = icxx"$(particle)->getEnergy();"
        SVector{4, Float64}(unsafe_load(p3, 1), unsafe_load(p3, 2), unsafe_load(p3, 3), e)
end

getType(particle) = icxx"$(particle)->getType();"


include("MCParticle.jl")
include("CaloHit.jl")
export CalHit, getP4, getPosition, CellIDDecoder,
    getEventNumber, getRunNumber, getDetectorName, getCollection, getCollectionNames, # LCEvent
    getType
#    getTypeName, # LCCollection
#    getEnergy, getParents, getDaughters, getPDG, getGeneratorStatus, getSimulatorStatus, isCreatedInSimulation, isBackScatter, vertexIsNotEndpointOfParent, #isDecayedInCalorimeter, hasLeftDetector, isStopped, isOverlay, getVertex, getTime, getEndpoint, getMomentum, getMomentumAtEndpoint, getMass, getCharge, # MCParticle
#    getCalorimeterHits, # Cluster
#    getClusters, getType, isCompound, getMass, getCharge, getReferencePoint, getParticleIDs, getParticleIDUsed, getGoodnessOfPID, getParticles, getClusters, getTracks, #getStartVertex, getEndVertex # ReconstructedParticle

#immutable CalHit
#       x::Cfloat
#       y::Cfloat
#       z::Cfloat
#       E::Cfloat
#end

#const MCPARTICLE = "MCParticle"
#const WRITE_NEW = 0
#const WRITE_APPEND = 1

# iteration over std vectors
#const StdVecs = Union{ClusterVec, CalorimeterHitVec, TrackVec, StringVec, MCParticleVec}

# uses Julia counting, 1..n
#start(it::StdVecs) = convert(UInt64, 1)
#next(it::StdVecs, i) = (it[i], i+1)
#done(it::StdVecs, i) = i > length(it)
#length(it::StdVecs) = size(it)
# 'at' uses C counting, 0..n-1
#getindex(it::StdVecs, i) = at(it, i-1)

#start(it::LCReader) = getNumberOfEvents(it)
#next(it::LCReader, state) = readNextEvent(it), state-1
#done(it::LCReader, state) = state < 1
#length(it::LCReader) = getNumberOfEvents(it)

#function iterate(f::Function, fn::AbstractString)
#    reader = createLCReader()
#    openFile(reader, fn)
#    try
#        for event in reader
#            f(event)
#        end
#    finally
#        closeFile(reader)
#        deleteLCReader(reader)
#    end
#end



# map from names stored in collection to actual types
#LCIOTypemap = Dict(
#    "CalorimeterHit" => CalorimeterHit,
#    "Cluster" => Cluster,
#       "LCGenericObject" => LCGenericObject,
#    "LCRelation" => LCRelation,
#       "MCParticle" => MCParticle,
#    "RawCalorimeterHit" => RawCalorimeterHit,
#    "ReconstructedParticle" => ReconstructedParticle,
#       "SimCalorimeterHit" => SimCalorimeterHit,
#       "SimTrackerHit" => SimTrackerHit,
#       "Track" => Track,
#       "TrackerHit" => TrackerHit,
#    "TrackerRawData" => TrackerRawData,
#    "Vertex" => Vertex,
#)

# This version of the iteration runs length() multiple times during the iteration
# if this becomes a speed problem, the length could be memoized, or iteration order could be inverted
#start(it::TypedCollection) = convert(UInt64, 1)
#done(it::TypedCollection, i) = i > length(it)
#next{T}(it::TypedCollection{T}, i) = it[i], i+1
#length(it::TypedCollection) = getNumberOfElements(it)
# getindex uses Julia counting, getElementAt uses C counting
#getindex(it::TypedCollection, i) = getElementAt(it, convert(UInt64, i-1))

#CellIDDecoder{T}(t::TypedCollection{T}) = CellIDDecoder{T}(coll(t))

#getTypeName{T}(coll::TypedCollection{T}) = "$T"

# to get the typed collection, one needs to read the typename
# then we can return the right type from the LCIOTypemap
#function getCollection(event, collectionName)
#       collection = getEventCollection(event, collectionName)
#       collectionType = getTypeName(collection)
#       return TypedCollection{LCIOTypemap[collectionType]}(collection)
#end

#type LCStdHepRdr
#    r::_LCStdHepRdrCpp
#    e::LCEventImpl
#end
#LCStdHepRdr(filename) = LCStdHepRdr(_LCStdHepRdrCpp(filename), LCEventImpl())

# the LCStdHepRdr implementation in C++ is not consistent with the LCIO reader
# this is me trying to make things a bit better
#start(it::LCStdHepRdr) = length(it.r)
#next(it::LCStdHepRdr, state) = readNextEvent(it.r), state-1
#done(it::LCStdHepRdr, state) = state < 1
#length(it::LCStdHepRdr) = getNumberOfEvents(it.r)

#function openStdhep(f::Function, fn::AbstractString)
#    reader = LCStdHepRdr(string(fn))
#    try
#       f(reader)
#    end
#end

# stdhep files only contain one specific type of collection: MCParticle
#function readNextEvent(r::LCStdHepRdr)
#    updateNextEvent(r.r, r.e)
#    getCollection(r.e, "MCParticle")
#end
#function getPosition(hit)
#    p3 = Array{Float64,1}(3)
#    valid = getPosition3(hit, p3)
#    return p3
#end

#function getMomentum(particle)
#    p3 = Array{Float64,1}(3)
#    valid = getMomentum3(particle, p3)
#    return p3
#end

#function getVertex(particle)
#    p3 = Array{Float64,1}(3)
#    valid = getVertex3(particle, p3)
#    return valid, p3
#end

#function getEndpoint(particle)
#    p3 = Array{Float64,1}(3)
#    valid = getEndpoint3(particle, p3)
#    return valid, p3
#end

#function getMomentumAtEndpoint(particle)
#    p3 = Array{Float64,1}(3)
#    valid = getMomentumAtEndpoint3(particle, p3)
#    return valid, p3
#end

# the navigator gets initialized with a collection
# it defers the actual work to the C++ implementation
#immutable LCRelationNavigator
#    relnav
#    fromType
#    toType
#    LCRelationNavigator(coll::TypedCollection) = _completNavigator(new(LCRelNav(coll.coll)))
#end
#function _completNavigator(nav)
#    nav.fromType = LCIOTypemap[nav.relnav.getFromType()]
#    nav.toType = LCIOTypemap[nav.relnav.getToType()]
#    nav
#end

# this ensures that the types are appropriately cast
#function getRelatedToObjects(nav::LCRelationNavigator, obj)
#    [CastOperator{nav.toType}.cast(x) for x in getRelatedToObjects(nav.relnav)]
#end

# this ensures that the types are appropriately cast
#function getRelatedFromObjects(nav::LCRelationNavigator, obj)
#    [CastOperator{nav.fromType}.cast(x) for x in getRelatedFromObjects(nav.relnav)]
#end

# should work for all particle types
#function getP4(x)
#    p3 = getMomentum(x)
#    E = getEnergy(x)
#    return (E, p3)
#end

#getP3(x) = getPosition(x)

# converters to keep older code working
#const CalHits = Union{SimCalorimeterHit, CalorimeterHit, RawCalorimeterHit}
#
#function CalHit(h::CalHits)
#    p = getPosition(h)
#    E = getEnergy(h)
#    return CalHit(p[1], p[2], p[3], E)
#end
#function convert(::Type(CalHit), h::CalHits)
#    p = getPosition(h)
#    E = getEnergy(h)
#    return CalHit(p[1], p[2], p[3], E)
#end

#function printParameters(p::LCParameters)
#    println("strings:")
#    for k in getStringKeys(p, StringVec())
#        println(k, "\t", getStringVal(p, k))
#    end
#    println("floats:")
#    for k in getFloatKeys(p, StringVec())
#        println(k, "\t", getFloatVal(p, k))
#    end
#    println("ints:")
#    for k in getIntKeys(p, StringVec())
#        println(k, "\t", getIntVal(p, k))
#    end
#end

end # module
