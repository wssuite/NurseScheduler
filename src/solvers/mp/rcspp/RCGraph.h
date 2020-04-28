//
// Created by antoine legrain on 2020-04-10.
//

#ifndef NURSESCHEDULER_RCGRAPH_H
#define NURSESCHEDULER_RCGRAPH_H

#include "tools/MyTools.h"
#include <boost/graph/adjacency_list.hpp>
#include "boost/config.hpp"
#include <boost/graph/r_c_shortest_paths.hpp>


enum LABEL {MAX_CONS_DAYS = 0, MIN_CONS_DAYS = 1, MAX_DAYS = 2, MIN_DAYS = 3, MAX_WEEKEND = 4};
// true if the dominance is done with a descending order (lower the better), false otherwise
static const std::vector<bool> labelsOrder = { true, false, true, false, true };
static const std::vector<std::string> labelName = {
    "MAX_CONS_DAYS", "MIN_CONS_DAYS", "MAX_DAYS     ", "MIN_DAYS     ", "MAX_WEEKEND  "};

// Different node types and their names
//
enum NodeType {
    SOURCE_NODE, PRINCIPAL_NETWORK, PRICE_LABEL_ENTRANCE,
    PRICE_LABEL, PRICE_LABEL_EXIT, SINK_NODE,
    NONE_NODE};
static const std::vector<std::string> nodeTypeName = {
    "SOURCE_NODE", "PPL_NETWORK", "PRI_LAB__IN",
    "PRICE_LABEL", "PRI_LAB_OUT", "SINK_NODE  ",
    "NONE       "};


// Different arc types and their names
//
enum ArcType{
    SOURCE_TO_PRINCIPAL, SHIFT_TO_NEWSHIFT, SHIFT_TO_SAMESHIFT, SHIFT_TO_ENDSEQUENCE,
    REPEATSHIFT, PRINCIPAL_TO_PRICE_LABEL, PRICE_LABEL_IN_TO_PRICE_LABEL, PRICE_LABEL_TO_PRICE_LABEL_OUT,
    PRICE_LABEL_OUT_TO_SINK, NONE_ARC};
static const std::vector<std::string> arcTypeName = {
    "SOURCE_TO_PPL  ", "SHIFT_TO_NEWSH ", "SHIFT_TO_SAMESH", "SHIFT_TO_ENDSEQ",
    "REPEATSHIFT    ", "PPL_TO_PRI_LAB ", "PLBIN_TO_PRILAB", "PRILAB_TO_PLOUT",
    "PLOUT  TO  SINK", "NONE           "};

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
                       std::vector<int> lbs={}, std::vector<int> ubs={}, bool hard_lbs=false, bool forbidden=false ) :
        num( n ), type( t ), lbs( lbs ), ubs( ubs ), hard_lbs_(hard_lbs), forbidden(forbidden) {}

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

    // hard lb
    bool hard_lbs_;

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
    Arc_Properties( int n = 0, int origin=-1, int destination=-1, ArcType ty = NONE_ARC, double c = 0, std::vector<int> consumptions={},
        int day=-1, std::vector<int> shifts={}, bool forbidden=false) :
        num( n ), origin(origin), destination(destination),
        type(ty), cost( c ), initialCost( c ), consumptions( consumptions ),
        day( day ), shifts( shifts ), forbidden(forbidden) {}

    // id
    //
    int num;

    // VERTICES
    int origin, destination;

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
typedef boost::graph_traits<Graph>::vertex_descriptor vertex;
typedef boost::graph_traits<Graph>::edge_descriptor edge;

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


    int first_day = -1, day = -1;
#ifdef DBG
    int pred_arc = -1;
    std::vector<int> shifts_;
#endif

    int label_value(int l) const {
      return label_values.at(l);
    }

    int size() const {
      return label_values.size();
    }

#ifdef DBG
    void print() const;
#endif
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
    dominance_spp(const std::vector<bool>& labelOrder, const std::vector<int>& labelsMinLevel):
      labelsOrder_(labelOrder), labelsMinLevel_(labelsMinLevel) {}
    inline bool operator()( const spp_res_cont& res_cont_1, const spp_res_cont& res_cont_2 ) const;
  private:
    std::vector<bool> labelsOrder_;
    // if label value is below (or above based on order) this level, label cannot be dominated
    // this is necessary as labels pricing could be even 0 for dominated label and label  could be reseted
    std::vector<int> labelsMinLevel_;
};


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

// comparator to order the processing of the nodes in boost rc spp
typedef ks_smart_pointer<boost::r_c_shortest_paths_label<Graph, spp_res_cont> > Spplabel;
struct SpplabelComparator {
    bool operator()(const Spplabel& splabel1, const Spplabel& splabel2) const;
};

struct RCSolution {
    RCSolution(int firstDay, const std::vector<int>& shifts, double c=0):
        firstDay(firstDay), shifts(shifts), cost(c) {};
    RCSolution(double c=0): firstDay(-1), cost(c) {};

    int firstDay;
    std::vector<int> shifts;
    double cost = 0;

    std::string toString(std::vector<int> shiftIDToShiftTypeID={}) const;
};

// modified boost::r_c_shortest_paths_dispatch function
template<class Graph,
    class VertexIndexMap,
    class Resource_Container,
    class Resource_Extension_Function,
    class Dominance_Function,
    class Label_Allocator,
    class Visitor,
    class Splabel_Comparator>
void r_c_shortest_paths_dispatch(const Graph& g,
                                 const VertexIndexMap& vertex_index_map,
                                 vertex s,
                                 std::vector<vertex> t,
                                 vector2D<edge>& opt_solutions_spp,
                                 std::vector<spp_res_cont>& pareto_opt_rcs_spp,
                                 const Resource_Container& rc,
                                 const Resource_Extension_Function& ref,
                                 const Dominance_Function& dominance,
                                 Label_Allocator /*la*/,
                                 Visitor vis,
                                 Splabel_Comparator);

//---------------------------------------------------------------------------
//
// C l a s s e s   G r a p h s
//
// Contains a interface for subgraphs and the main rcspp
//
//---------------------------------------------------------------------------
class SubGraph {
  public:
    SubGraph() {};
    virtual ~SubGraph() {};

    virtual int entrance(int day=-1) const=0;
    virtual int exit(int day=-1) const=0;

    // link two sub graphs together:
    // 1. create an arc from the exit of inSubGraph to the current entrance
    // 2. the entrance method should now returns the entrance of the  inSubGraph
    virtual void linkInSubGraph(SubGraph& inSubGraph, int day=-1)=0;

    // link two sub graphs together:
    // 1. create an arc from the current exit to the entrance of outSubGraph
    // 2. the exit method should now returns the exit of the outSubGraph
    virtual void linkOutSubGraph(SubGraph& outSubGraph, int day=-1)=0;
};

class RCGraph {
  public:
    RCGraph(int nDays=0);
    virtual ~RCGraph();

    std::vector<RCSolution> solve(int nLabels, double maxReducedCostBound,
                                  const std::vector<int>& labelsMinLevel,
                                  std::vector<vertex> sinks={},
                                  std::function<void (spp_res_cont&)> post_process_rc=nullptr);

    ////////////////////
    //  NODES
    ////////////////////

    Graph& g() { return g_; }

    // Basic function for adding a node
    int addSingleNode(NodeType type, std::vector<int> lbs, std::vector<int> ubs, bool hard_lbs=false);
    void setSource(int v) { source_ = v; }
    int source() const { return source_; }
    void addSink(int v) { sinks_.push_back(v); }
    int sink(int k=0) const { return sinks_.at(k); }
    int lastSink() const { return sinks_.back(); }
    const std::vector<vertex>& sinks() const { return sinks_; }

    // Get info from the node ID
    inline int nodesSize() const { return nNodes_; }
    inline const Vertex_Properties & node(int v) const {
      if(v == -1) Tools::throwError("Cannot retrieve a node for index -1;");
      return get( boost::vertex_bundle, g_ )[v]; }
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
      if(a == -1) Tools::throwError("Cannot retrieve an arc for index -1;");
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
    void printGraph(int nLabel=-1, int nShifts=1) const;
    std::string printNode(int v, int nLabel=-1) const ;
    void printAllNodes(int nLabel=-1) const ;
    std::string printArc(int a, int nLabel=-1, int nShifts=1) const ;
    void printAllArcs(int nLabel=-1, int nShifts=1) const;
    std::string shortNameNode(int v) const;
    std::string printSummaryOfGraph() const;
    void printPath(std::ostream& out, std::vector<edge> path, spp_res_cont resource) const;

  protected:
    // THE GRAPH
    Graph g_;
    int nDays_;
    int nNodes_;										// Total number of nodes in the rcspp
    // Source
    vertex source_;
    // Sink Node
    std::vector<vertex> sinks_;

    int nArcs_;											// Total number of arcs in the rcspp
    std::vector<edge> arcsDescriptors_;

    std::set<int> forbiddenNodes_;
    std::set<int> forbiddenArcs_;

    bool processPath(std::vector<edge>& path, spp_res_cont& rc,
        std::function<void (spp_res_cont&)> post_process_rc) const;

    RCSolution solution(const std::vector<edge>& path,
        const spp_res_cont& resource) const;
};




#endif //NURSESCHEDULER_RCGRAPH_H
