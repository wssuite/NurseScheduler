//
// Created by antoine legrain on 2020-04-10.
//

#ifndef NURSESCHEDULER_RCGRAPH_H
#define NURSESCHEDULER_RCGRAPH_H

//#include "tools/MyTools.h"
#include <boost/graph/adjacency_list.hpp>
#include "boost/config.hpp"
#include <boost/graph/r_c_shortest_paths.hpp>


enum LABEL {MAX_CONS_DAYS = 0, MIN_CONS_DAYS = 1};
// true if the dominance is done with a descending order (lower the better), false otherwise
static const std::vector<bool> labelsOrder = { true, false };
static const std::vector<std::string> labelName = {
    "MAX_CONS_DAYS", "MIN_CONS_DAYS", "MAX_DAYS     ", "MIN_DAYS     "};

// Different node types and their names
//
enum NodeType {
    SOURCE_NODE, PRINCIPAL_NETWORK, ROTATION_LENGTH_ENTRANCE,
    ROTATION_LENGTH, ROTATION_LENGTH_EXIT, SINK_NODE,
    NONE_NODE};
static const std::vector<std::string> nodeTypeName = {
    "SOURCE_NODE", "PPL_NETWORK", "ROTSIZE__IN",
    "ROTSIZE    ", "ROTSIZE_OUT", "SINK_NODE  ",
    "NONE       "};


// Different arc types and their names
//
enum ArcType{
    SOURCE_TO_PRINCIPAL, SHIFT_TO_NEWSHIFT, SHIFT_TO_SAMESHIFT, SHIFT_TO_ENDSEQUENCE,
    REPEATSHIFT, PRINCIPAL_TO_ROTSIZE, ROTSIZEIN_TO_ROTSIZE, ROTSIZE_TO_ROTSIZEOUT,
    ROTSIZEOUT_TO_SINK, NONE_ARC};
static const std::vector<std::string> arcTypeName = {
    "SOURCE_TO_PPL  ", "SHIFT_TO_NEWSH ", "SHIFT_TO_SAMESH", "SHIFT_TO_ENDSEQ",
    "REPEATSHIFT    ", "PPL_TO_ROTSIZE ", "ROTSZIN_TO_RTSZ", "ROTSIZE_TO_ROUT",
    "ROT OUT TO SINK", "NONE           "};

//////////////////////////////////////////////////////////////////////////
//
// The following structures are the properties of the nodes and arcs:
//   For the nodes: id, earliest arrival time, and latest arrival time
//   For the arcs : id, cost and time (time is the resource here)
//   For the rcspp: usual stuff + nodes and and special properties
//
//////////////////////////////////////////////////////////////////////////

// Nodes specific properties for RC
//
struct Vertex_Properties{

    // Constructor
    //
    Vertex_Properties( int n = 0, NodeType t = NONE_NODE,
                       std::vector<int> lbs={}, std::vector<int> ubs={}, bool forbidden=false ) :
        num( n ), type( t ), lbs( lbs ), ubs( ubs ), forbidden(forbidden) {}

    // id
    //
    int num;

    // type
    //
    NodeType type;

    // label lower bounds
    //
    std::vector<int> lbs;

    // label upper bounds
    //
    std::vector<int> ubs;

    // forbidden
    //
    bool forbidden;

    int lb(int l) const {
      return lbs.at(l);
    }

    int ub(int l) const {
      return ubs.at(l);
    }

    int size() const {
      return lbs.size();
    }
};

// Arcs specific properties for RC
//
struct Arc_Properties{

    // Constructor
    //
    Arc_Properties( int n = 0, ArcType ty = NONE_ARC, double c = 0, std::vector<int> consumptions={},
        int day=-1, std::vector<int> shifts={}, bool forbidden=false) :
        num( n ), type(ty), cost( c ), initialCost( c ), consumptions( consumptions ),
        day( day ), shifts( shifts ), forbidden(forbidden) {}

    // id
    //
    int num;

    // type
    //
    ArcType type;

    // traversal cost
    //
    double cost;
    double initialCost;

    // label consumption
    std::vector<int> consumptions;

    int day; // day
    std::vector<int> shifts; // shifts id

    bool forbidden;

    int consumption(int l) const {
      return consumptions.at(l);
    }

    int size() const {
      return consumptions.size();
    }
};

// Graph with RC generic structure
//
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, Vertex_Properties, Arc_Properties> Graph;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// The following functions / data structures are used for the time resource:
//   Resource container ("resource" = cost + time)
//   Comparison override --> in SubProblem.cpp
//   Cost extension function
//   Dominance function
//
/////////////////////////////////////////////////////////////////////////////

// data structures for shortest path problem with labels
// ResourceContainer model

struct spp_res_cont{

    // Constructor
    //
    spp_res_cont( double c, const std::vector<int>& label_values ) :
        cost( c ), label_values( label_values ) {}

    // Assign
    //
    spp_res_cont& operator=( const spp_res_cont& other ){
      if( this == &other ) return *this;
      this->~spp_res_cont();
      new( this ) spp_res_cont( other );
      return *this;
    }

    // Current cost
    //
    double cost;

    // Current labels
    //
    std::vector<int> label_values;

    int label_value(int l) const {
      return label_values.at(l);
    }

    int size() const {
      return label_values.size();
    }
};

// Resources extension model (arc has cost + label consumptions)
class ref_spp{
  public:
    inline bool operator()( const Graph& g,
                            spp_res_cont& new_cont,
                            const spp_res_cont& old_cont,
                            boost::graph_traits<Graph>::edge_descriptor ed ) const;
};

// Dominance function model
class dominance_spp{
  public:
    inline bool operator()( const spp_res_cont& res_cont_1, const spp_res_cont& res_cont_2 ) const;
};

//----------------------------------------------------------------
//
// Shortest path function with several sinks
// (modified from boost so that we can give several sink nodes)
//
//----------------------------------------------------------------

// r_c_shortest_paths_dispatch function (body/implementation)
template<class Graph,
    class VertexIndexMap,
    class EdgeIndexMap,
    class Resource_Container,
    class Resource_Extension_Function,
    class Dominance_Function,
    class Label_Allocator,
    class Visitor>
void r_c_shortest_paths_dispatch_several_sinks
    ( const Graph& g,
      const VertexIndexMap& vertex_index_map,
      const EdgeIndexMap& /*edge_index_map*/,
      typename boost::graph_traits<Graph>::vertex_descriptor s,
      std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> t,
        // each inner vector corresponds to a pareto-optimal path
      std::vector
      <std::vector
          <typename boost::graph_traits
              <Graph>::edge_descriptor> >& pareto_optimal_solutions,
      std::vector
      <Resource_Container>& pareto_optimal_resource_containers,
      bool b_all_pareto_optimal_solutions,
        // to initialize the first label/resource container
        // and to carry the type information
      const Resource_Container& rc,
      const Resource_Extension_Function& ref,
      const Dominance_Function& dominance,
        // to specify the memory management strategy for the labels
      Label_Allocator /*la*/,
      Visitor vis );

/////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
//
// C l a s s   k s _ s m a r t _ p o i n t e r
//
// Needed for the shortest path algorithms.
// From: Kuhlin, S.; Schader, M. (1999): Die C++-Standardbibliothek, Springer, Berlin, p. 333 f.
//
//---------------------------------------------------------------------------
template<class T>
class ks_smart_pointer
{
  public:
    ks_smart_pointer( T* ptt = 0 ) : pt( ptt ) {}
    ks_smart_pointer( const ks_smart_pointer& other ) : pt( other.pt ) {}
    ks_smart_pointer& operator=( const ks_smart_pointer& other )
    { pt = other.pt; return *this; }
    ~ks_smart_pointer() {}
    T& operator*() const { return *pt; }
    T* operator->() const { return pt; }
    T* get() const { return pt; }
    operator T*() const { return pt; }
    friend bool operator==( const ks_smart_pointer& t,
                            const ks_smart_pointer& u )
    { return *t.pt == *u.pt; }
    friend bool operator!=( const ks_smart_pointer& t,
                            const ks_smart_pointer& u )
    { return *t.pt != *u.pt; }
    friend bool operator<( const ks_smart_pointer& t,
                           const ks_smart_pointer& u )
    { return *t.pt < *u.pt; }
    friend bool operator>( const ks_smart_pointer& t,
                           const ks_smart_pointer& u )
    { return *t.pt > *u.pt; }
    friend bool operator<=( const ks_smart_pointer& t,
                            const ks_smart_pointer& u )
    { return *t.pt <= *u.pt; }
    friend bool operator>=( const ks_smart_pointer& t,
                            const ks_smart_pointer& u )
    { return *t.pt >= *u.pt; }
  private:
    T* pt;
}; // ks_smart_pointer


struct RCSolution {
    RCSolution(int firstDay, const std::vector<int>& shifts, double c=0):
        firstDay(firstDay), shifts(shifts), cost(c) {};
    RCSolution(double c=0): firstDay(-1), cost(c) {};

    int firstDay;
    std::vector<int> shifts;
    double cost = 0;

    std::string toString(std::vector<int> shiftIDToShiftTypeID={}) const;
};

//---------------------------------------------------------------------------
//
// C l a s s e s   G r a p h s
//
// Contains the main rcspp and some reccurent subgraphs
//
//---------------------------------------------------------------------------

class RCGraph {
  public:
    RCGraph(int nDays=0);
    virtual ~RCGraph();

    std::vector<RCSolution> solve(int nLabels, double maxReducedCostBound,
        std::vector<boost::graph_traits<Graph>::vertex_descriptor> sinks={});

    RCSolution solution(
        const std::vector< boost::graph_traits<Graph>::edge_descriptor >& path,
        const spp_res_cont& resource);

    ////////////////////
    //  NODES
    ////////////////////

    // Basic function for adding a node
    int addSingleNode(NodeType type, std::vector<int> lbs, std::vector<int> ubs);
    void setSource(int v) { source_ = v; }
    int source() const { return source_; }
    void addSink(int v) { sinks_.push_back(v); }
    int sink(int k=0) const { return sinks_.at(k); }
    int lastSink() const { return sinks_.back(); }
    const std::vector<boost::graph_traits<Graph>::vertex_descriptor>& sinks() const { return sinks_; }

    // Get info from the node ID
    inline int nodesSize() const { return nNodes_; }
    inline const Vertex_Properties & node(int v) const { return get( boost::vertex_bundle, g_ )[v]; }
    inline NodeType nodeType(int v) const {return get( &Vertex_Properties::type, g_)[v];}
    inline const std::vector<int> & nodeLBs(int v) const {return get( &Vertex_Properties::lbs, g_)[v];}
    inline const std::vector<int> & nodeUBs(int v) const {return get( &Vertex_Properties::ubs, g_)[v];}
    inline bool nodeForbidden(int v) const {return forbiddenNodes_.find(v) != forbiddenNodes_.end();}

    inline void updateUBs(int v, const std::vector<int>& ubs){boost::put( &Vertex_Properties::ubs, g_, v, ubs);}
    inline void forbidNode(int v) {
      boost::put( &Vertex_Properties::forbidden, g_, v, true);
      forbiddenNodes_.insert(v);
    }
    inline void authorizeNode(int v) {
      boost::put( &Vertex_Properties::forbidden, g_, v, false);
      forbiddenNodes_.erase(v);
    }

    ////////////////////
    //  ARCS
    ////////////////////

    // Basic function for adding an arc
    int addSingleArc(int origin, int destination, double baseCost, std::vector<int> consumptions,
        ArcType type, int day = -1, std::vector<int> shifts = {});

    // Get info with the arc ID
    inline int arcsSize() const { return nArcs_; }
    inline const Arc_Properties & arc(int a) const {
      return get( boost::edge_bundle, g_ )[arcsDescriptors_[a]];
    }
    inline ArcType arcType(int a) const {return get( &Arc_Properties::type, g_, arcsDescriptors_[a]);}
    inline int arcOrigin(int a) const {return boost::source(arcsDescriptors_[a], g_);}
    inline int arcDestination(int a) const {return target(arcsDescriptors_[a], g_);}
    inline const std::vector<int> & arcConsumptions(int a) const {
      return get( &Arc_Properties::consumptions, g_, arcsDescriptors_[a]);
    }
    inline double arcCost(int a) const {return get( &Arc_Properties::cost, g_, arcsDescriptors_[a]);}
    inline double arcInitialCost(int a) const {return get( &Arc_Properties::initialCost, g_, arcsDescriptors_[a]);}
    inline const std::vector<int>& arcShifts(int a) const {return get( &Arc_Properties::shifts, g_, arcsDescriptors_[a]);}
    inline int arcDay(int a) const {return get( &Arc_Properties::day, g_, arcsDescriptors_[a]);}
    inline bool arcForbidden(int a) const {return get( &Arc_Properties::forbidden, g_, arcsDescriptors_[a]);}

    inline void updateConsumptions(int a, const std::vector<int>& consumptions){
      boost::put( &Arc_Properties::consumptions, g_, arcsDescriptors_[a], consumptions );
    }
    inline void updateShifts(int a, const std::vector<int>& shifts){
      boost::put( &Arc_Properties::shifts, g_, arcsDescriptors_[a], shifts );
    }
    inline void updateCost(int a, double cost){
      boost::put( &Arc_Properties::cost, g_, arcsDescriptors_[a], cost );
    }
    inline void forbidArc(int a) {
      boost::put( &Arc_Properties::forbidden, g_, arcsDescriptors_[a], true);
      forbiddenArcs_.insert(a);
    }
    inline void authorizeArc(int a) {
      boost::put( &Arc_Properties::forbidden, g_, arcsDescriptors_[a], false);
      forbiddenArcs_.erase(a);
    }

    void resetAuthorizations();

    // Print functions.
    //
    void printGraph() const;
    std::string printNode(int v) const ;
    void printAllNodes() const ;
    std::string printArc(int a) const ;
    void printAllArcs() const;
    std::string shortNameNode(int v) const;
    std::string printSummaryOfGraph() const;
    void printPath(std::vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_res_cont ressource) const;

    // Test function for Shortest Path Problem with Resource Constraint
    //
    void testGraph_spprc();

  protected:
    // THE GRAPH
    Graph g_;
    int nDays_;
    int nNodes_;										// Total number of nodes in the rcspp
    // Source
    typename boost::graph_traits<Graph>::vertex_descriptor source_;
    // Sink Node
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> sinks_;

    int nArcs_;											// Total number of arcs in the rcspp
    std::vector< boost::graph_traits< Graph>::edge_descriptor > arcsDescriptors_;

    std::set<int> forbiddenNodes_;
    std::set<int> forbiddenArcs_;
};


#endif //NURSESCHEDULER_RCGRAPH_H
