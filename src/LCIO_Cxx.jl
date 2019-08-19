__precompile__(false)
module LCIO_Cxx
using Cxx
using StaticArrays
using Libdl

import Base: getindex, iterate, length, convert, +

const lcio_path="/home/jstrube/.julia/dev/LCIO/deps/usr/"
addHeaderDir(lcio_path*"include", kind=C_System)
Libdl.dlopen(lcio_path*"lib/liblcio.so", Libdl.RTLD_GLOBAL)
cxx"""
#include "lcio.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "UTIL/CellIDDecoder.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/Track.h"
#include "EVENT/TrackerRawData.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"
#include "EVENT/LCRelation.h"
#include "EVENT/LCGenericObject.h"
#include "EVENT/Cluster.h"
#include "EVENT/CalorimeterHit.h"
#include <vector>
#include <iostream>
#include <string>
"""

function __init__()
end

const ClusterVec = cxxt"EVENT::ClusterVec"
const CalorimeterHitVec = cxxt"EVENT::CalorimeterHitVec"
const TrackVec = cxxt"EVENT::TrackVec"
const StringVec = cxxt"EVENT::StringVec"
const MCParticleVec = cxxt"EVENT::MCParticleVec"

# iteration over std vectors
const StdVecs = Union{ClusterVec, CalorimeterHitVec, TrackVec, StringVec, MCParticleVec}

# uses Julia counting, 1..n
iterate(it::StdVecs) = length(it) > 0 ? (it[1], 2) : nothing
iterate(it::StdVecs, i) = i <= length(it) ? (it[i], i+1) : nothing
length(it::StdVecs) = icxx"$(it)->size;"
# 'at' uses C counting, 0..n-1
# FIXME is the cast necessary?
getindex(it::StdVecs, i) = at(it, convert(UInt64, i-1))
eltype(::Type{ClusterVec}) = Cluster
eltype(::Type{CalorimeterHitVec}) = CalorimeterHit
eltype(::Type{TrackVec}) = Track
eltype(::Type{StringVec}) = String
eltype(::Type{MCParticleVec}) = MCParticle

const LCReader = cxxt"IO::LCReader*"
function iterate(it::LCReader)
    event = icxx"$(it)->readNextEvent();"
    event == C_NULL && return nothing
    return (event, nothing)
end
iterate(it::LCReader, state) = iterate(it)
length(it::LCReader) = icxx"$(it)->getNumberOfEvents();"
eltype(::Type{LCReader}) = LCEvent
    
function open(f::Function, fn::AbstractString)
    reader = icxx"IOIMPL::LCFactory::getInstance()->createLCReader();"
    reader == C_NULL && return nothing
    try
        icxx"""$(reader)->open($(fn));"""
        f(reader)
    finally
        icxx"""$(reader)->close();"""
        icxx"""delete $(reader);"""
    end
end
    
    
# The iterators use the events as state
# the lcio reader returns NULL at the end of the file
# Because the LCIO getNext function could re-use the memory for the next event,
# the next function should not hold the current and next event at the same time.
# Returning current and reading nextEvent as the new state causes memory corruption
const LCEvent = cxxt"EVENT::LCEvent*"
# uses Julia counting, 1..n
iterate(it::StdVecs) = length(it) > 0 ? (it[1], 2) : nothing
iterate(it::StdVecs, i) = i <= length(it) ? (it[i], i+1) : nothing
length(it::StdVecs) = size(it)
# 'at' uses C counting, 0..n-1
# FIXME is the cast necessary?
getindex(it::StdVecs, i) = at(it, convert(UInt64, i-1))

function iterate(it::LCEvent)
    event = icxx"$(reader)->readNextEvent();"
    event != C_NULL ? (event, nothing) : nothing
end
iterate(it::LCEvent, state) = iterate(it)
length(it::LCEvent) = icxx"$(reader)->getNumberOfEvents();"

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
# open file with reader, returns iterator

# We would like to have a typed collection, but what getCollection returns is unfortunately untyped
# The type is established by reading its name from the collection and mapping it in the LCIOTypemap
struct LCCollection{T}
        coll::cxxt"EVENT::LCCollection*"
end

const CalorimeterHit = cxxt"EVENT::CalorimeterHit*"
const Cluster = cxxt"EVENT::Cluster*"
const LCGenericObject = cxxt"EVENT::LCGenericObject*"
const LCRelation = cxxt"EVENT::LCRelation*"
const MCParticle = cxxt"EVENT::MCParticle*"
const ReconstructedParticle = cxxt"EVENT::ReconstructedParticle*"
const SimCalorimeterHit = cxxt"EVENT::SimCalorimeterHit*"
const SimTrackerHit = cxxt"EVENT::SimTrackerHit*"
const Track = cxxt"EVENT::Track*"
const TrackerHit = cxxt"EVENT::TrackerHit*"
const TrackerRawData = cxxt"EVENT::TrackerRawData*"
const Vertex = cxxt"EVENT::Vertex*"

# map from names stored in collection to actual types
const LCIOTypemap = Dict(
    "CalorimeterHit" => CalorimeterHit,
    "Cluster" => Cluster,
    "LCGenericObject" => LCGenericObject,
    "LCRelation" => LCRelation,
    "MCParticle" => MCParticle,
    "ReconstructedParticle" => ReconstructedParticle,
    "SimCalorimeterHit" => SimCalorimeterHit,
    "SimTrackerHit" => SimTrackerHit,
    "Track" => Track,
    "TrackerHit" => TrackerHit,
    "TrackerRawData" => TrackerRawData,
    "Vertex" => Vertex,
)

# CellIDDecoder{T}(t::LCCollection{T}) = icxx"UTIL::CellIDDecoder<$(T)>($(t.coll));"


decode(iddecoder) = icxx"$(iddecoder)->decode();"



iterate(it::LCCollection{T}) where {T} = length(it) > 0 ? (it[1], 2) : nothing
iterate(it::LCCollection{T}, i) where {T} = i <= length(it) ? (it[i], i+1) : nothing
length(it::LCCollection) = icxx"$(it.coll)->getNumberOfElements();"
getindex(it::LCCollection{T}, i) where {T} = icxx"static_cast<$(T)>($(it.coll)->getElementAt($(i)-1));"

function getCollection(event, collectionName)
    collection = icxx"""$(event)->getCollection($(pointer(collectionName)));"""
    collectionType = icxx"$(collection)->getTypeName();"
    return LCCollection{LCIOTypemap[String(collectionType)]}(collection)
end

getTypeName(collection::LCCollection) = String( icxx"$(collection.coll)->getTypeName();" )
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

getDetectorName(event) = String(icxx"$(event)->getDetectorName();")

include("MCParticle.jl")
include("CaloHit.jl")
export CalHit, getP4, getPosition, 
    getEventNumber, getRunNumber, getDetectorName, getCollection, getCollectionNames, # LCEvent
    getType, getEnergy, getMomentum, getPDG, getParents, getTypeName,
    decode
#    getTypeName, # LCCollection
#    getEnergy, getParents, getDaughters, getPDG, getGeneratorStatus, getSimulatorStatus, isCreatedInSimulation, isBackScatter, vertexIsNotEndpointOfParent, #isDecayedInCalorimeter, hasLeftDetector, isStopped, isOverlay, getVertex, getTime, getEndpoint, getMomentum, getMomentumAtEndpoint, getMass, getCharge, # MCParticle
#    getCalorimeterHits, # Cluster
#    getClusters, getType, isCompound, getMass, getCharge, getReferencePoint, getParticleIDs, getParticleIDUsed, getGoodnessOfPID, getParticles, getClusters, getTracks, #getStartVertex, getEndVertex # ReconstructedParticle


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


# This version of the iteration runs length() multiple times during the iteration
# if this becomes a speed problem, the length could be memoized, or iteration order could be inverted
#start(it::TypedCollection) = convert(UInt64, 1)
#done(it::TypedCollection, i) = i > length(it)
#next{T}(it::TypedCollection{T}, i) = it[i], i+1
#length(it::TypedCollection) = getNumberOfElements(it)
# getindex uses Julia counting, getElementAt uses C counting
#getindex(it::TypedCollection, i) = getElementAt(it, convert(UInt64, i-1))

#getTypeName{T}(coll::TypedCollection{T}) = "$T"

# to get the typed collection, one needs to read the typename
# then we can return the right type from the LCIOTypemap
#function getCollection(event, collectionName)
#       collection = getEventCollection(event, collectionName)
#       collectionType = getTypeName(collection)
#       return TypedCollection{LCIOTypemap[collectionType]}(collection)
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
