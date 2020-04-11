//
// Created by antoine legrain on 2020-04-10.
//

#ifndef NURSESCHEDULER_RCGRAPH_H
#define NURSESCHEDULER_RCGRAPH_H

#include "tools/MyTools.h"
#include <boost/graph/adjacency_list.hpp>
#include "boost/config.hpp"
#include <boost/graph/r_c_shortest_paths.hpp>


enum LABEL {MAX_CONS_DAYS = 0, MIN_CONS_DAYS = 1};
// true if the dominance is done with a descending order (lower the better), false otherwise
static const std::vector<bool> labelsOrder = { true, false };
static const std::vector<string> labelName = {
    "MAX_CONS_DAYS", "MIN_CONS_DAYS", "MAX_DAYS     ", "MIN_DAYS     "};

// Different node types and their names
//
enum NodeType {
    SOURCE_NODE, PRINCIPAL_NETWORK, ROTATION_LENGTH_ENTRANCE,
    ROTATION_LENGTH, ROTATION_LENGTH_EXIT, SINK_NODE,
    NONE_NODE};
static const std::vector<string> nodeTypeName = {
    "SOURCE_NODE", "PPL_NETWORK", "ROTSIZE__IN",
    "ROTSIZE    ", "ROTSIZE_OUT", "SINK_NODE  ",
    "NONE       "};


// Different arc types and their names
//
enum ArcType{
    SOURCE_TO_PRINCIPAL, SHIFT_TO_NEWSHIFT, SHIFT_TO_SAMESHIFT, SHIFT_TO_ENDSEQUENCE,
    REPEATSHIFT, PRINCIPAL_TO_ROTSIZE, ROTSIZEIN_TO_ROTSIZE, ROTSIZE_TO_ROTSIZEOUT,
    ROTSIZEOUT_TO_SINK, NONE_ARC};
static const std::vector<string> arcTypeName = {
    "SOURCE_TO_PPL  ", "SHIFT_TO_NEWSH ", "SHIFT_TO_SAMESH", "SHIFT_TO_ENDSEQ",
    "REPEATSHIFT    ", "PPL_TO_ROTSIZE ", "ROTSZIN_TO_RTSZ", "ROTSIZE_TO_ROUT",
    "ROT OUT TO SINK", "NONE           "};

//////////////////////////////////////////////////////////////////////////
//
// The following structures are the properties of the nodes and arcs:
//   For the nodes: id, earliest arrival time, and latest arrival time
//   For the arcs : id, cost and time (time is the resource here)
//   For the graph: usual stuff + nodes and and special properties
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
                            boost::graph_traits<Graph>::edge_descriptor ed ) const{
      const Arc_Properties& arc_prop = get( boost::edge_bundle, g )[ed];
      if(arc_prop.forbidden) return false;
      const Vertex_Properties& vert_prop = get( boost::vertex_bundle, g )[target( ed, g )];
      if(vert_prop.forbidden) return false;
      new_cont.cost = old_cont.cost + arc_prop.cost;

      for(int l=0; l<old_cont.size(); ++l) {
        int lv = std::max(vert_prop.lb(l), old_cont.label_value(l) + arc_prop.consumption(l));
        if(lv > vert_prop.ub(l))
          return  false;
        new_cont.label_values[l] = lv;
      }

      return true;
    }
};

// Dominance function model
class dominance_spp{
  public:
    inline bool operator()( const spp_res_cont& res_cont_1, const spp_res_cont& res_cont_2 ) const	{
      // must be "<=" here!!!
      // must NOT be "<"!!!
      if(res_cont_1.cost > res_cont_2.cost) return false;
      for(int l=0; l<res_cont_1.size(); ++l)
        if(labelsOrder.at(l)) { // dominance done with descending order (lower the better)
          if(res_cont_1.label_value(l) > res_cont_2.label_value(l)) return false;
        } else if(res_cont_1.label_value(l) < res_cont_2.label_value(l)) return false;
      return true;
      // this is not a contradiction to the documentation
      // the documentation says:
      // "A label $l_1$ dominates a label $l_2$ if and only if both are resident
      // at the same vertex, and if, for each resource, the resource consumption
      // of $l_1$ is less than or equal to the resource consumption of $l_2$,
      // and if there is at least one resource where $l_1$ has a lower resource
      // consumption than $l_2$."
      // one can think of a new label with a resource consumption equal to that
      // of an old label as being dominated by that old label, because the new
      // one will have a higher number and is created at a later point in time,
      // so one can implicitly use the number or the creation time as a resource
      // for tie-breaking
    }
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
    RCSolution(double c=0): firstDay(-1), cost(c) {};
    int firstDay;
    std::vector<int> shifts;
    double cost = 0;
};

//---------------------------------------------------------------------------
//
// C l a s s e s   G r a p h s
//
// Contains the main graph and some reccurent subgraphs
//
//---------------------------------------------------------------------------

class RCGraph {
  public:
    RCGraph(int nDays=0);
    virtual ~RCGraph();

    std::vector<RCSolution> solve(int nLabels, double maxReducedCostBound);

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
    inline ArcType arcType(int a) const {return allArcsTypes_[a];}
    inline int arcOrigin(int a) const {return boost::source(arcsDescriptors_[a], g_);}
    inline int arcDestination(int a) const {return target(arcsDescriptors_[a], g_);}
    inline const std::vector<int> & arcConsumptions(int a) const {
      return get( &Arc_Properties::consumptions, g_, arcsDescriptors_[a]);
    }
    inline double arcCost(int a) const {return get( &Arc_Properties::cost, g_, arcsDescriptors_[a]);}
    inline double arcInitialCost(int a) const {return get( &Arc_Properties::initialCost, g_, arcsDescriptors_[a]);}
    inline const std::vector<int>& arcShifts(int a) const {return get( &Arc_Properties::shifts, g_, arcsDescriptors_[a]);}
    inline int arcDay(int a) const {return get( &Arc_Properties::day, g_, arcsDescriptors_[a]);}
    inline bool arcForbidden(int a) const {return forbiddenArcs_.find(a) != forbiddenArcs_.end();}

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
    void printGraph();
    string printNode(int v);
    void printAllNodes();
    string printArc(int a);
    void printAllArcs();
    string shortNameNode(int v);
    string printSummaryOfGraph();
    virtual void printPath(std::vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_res_cont ressource);

    // Test function for Shortest Path Problem with Resource Constraint
    //
    void testGraph_spprc();

  protected:
    // THE GRAPH
    Graph g_;
    int nDays_;
    int nNodes_;										// Total number of nodes in the graph
    // Source
    typename boost::graph_traits<Graph>::vertex_descriptor source_;
    // Sink Node
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> sinks_;

    int nArcs_;											// Total number of arcs in the graph
    std::vector<ArcType> allArcsTypes_;						// Vector of their types
    std::vector< boost::graph_traits< Graph>::edge_descriptor > arcsDescriptors_;

    std::set<int> forbiddenNodes_;
    std::set<int> forbiddenArcs_;
};

//class PrincipalSubGraph {
//  public:
//    PrincipalSubGraph();
//    virtual ~PrincipalSubGraph();
//
//  protected:
//    // Nodes of the PRINCIPAL_NETWORK subnetwork
//    vector3D principalNetworkNodes_;					// For each SHIFT, DAY, and # of CONSECUTIVE, the corresponding node id
//    std::vector<int> maxvalConsByShift_;						// For each shift, number of levels that the subnetwork contains
//    std::map<int,int> principalToShift_;						// For each node of the principal network, maps it ID to the shift it represents
//    std::map<int,int> principalToDay_;						// For each node of the principal network, maps it ID to the day it represents
//    std::map<int,int> principalToCons_;						// For each node of the principal network, maps it ID to the number of consecutive shifts it represents
//
//    vector4D arcsShiftToNewShift_;		// Index: (shiftType1, shiftType2, day1, shift)
//    vector4D arcsShiftToSameShift_;		// Index: (shiftType, day, nCons, shift) of origin
//    vector3D arcsShiftToEndsequence_;	// Index: (shiftType, day, nCons) of origin
//    vector3D arcsRepeatShift_;		// Index: (shiftType, day, shift) of origin
//};
//
//class RotationLengthSubGraph {
//  public:
//    RotationLengthSubGraph();
//    virtual ~RotationLengthSubGraph();
//
//  protected:
//    // Nodes of the ROTATION_LENGTH subnetwork
//    std::vector<int> rotationLengthEntrance_;				// For each day, entrance node to the ROTATION_LENGTH subnetwork
//    std::vector<std::map<int,int> > rotationLengthNodes_;			// For each day, maps the length of the rotation to the corresponding check node
//    std::map<int,int> rotationLengthNodesLAT_;				// For each rotation length node, the corresponding EAT
//
//    std::vector<std::map<int,int> > arcsRotsizeinToRotsizeDay_;	// Index: (day,size) of the rotation [destination]
//    std::vector<std::map<int,int> > arcsRotsizeToRotsizeoutDay_;	// Index: (day,size) of the rotation [origin]
//};


#endif //NURSESCHEDULER_RCGRAPH_H
