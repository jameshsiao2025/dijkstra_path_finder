// Digraph.hpp
//
// ICS 46 Winter 2020
// Project #5: Rock and Roll Stops the Traffic
//
// This header file declares a class template called Digraph, which is
// intended to implement a generic directed graph.  The implementation
// uses the adjacency lists technique, so each vertex stores a linked
// list of its outgoing edges.
//
// Along with the Digraph class template is a class DigraphException
// and a couple of utility structs that aren't generally useful outside
// of this header file.
//
// In general, directed graphs are all the same, except in the sense
// that they store different kinds of information about each vertex and
// about each edge; these two types are the type parameters to the
// Digraph class template.

#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <exception>
#include <functional>
#include <list>
#include <map>
#include <utility>
#include <vector>

//#include <iostream>
#include <queue>



// DigraphExceptions are thrown from some of the member functions in the
// Digraph class template, so that exception is declared here, so it
// will be available to any code that includes this header file.

class DigraphException : public std::runtime_error
{
public:
    DigraphException(const std::string& reason);
};


inline DigraphException::DigraphException(const std::string& reason)
    : std::runtime_error{reason}
{
}



// A DigraphEdge lists a "from vertex" (the number of the vertex from which
// the edge points), a "to vertex" (the number of the vertex to which the
// edge points), and an EdgeInfo object.  Because different kinds of Digraphs
// store different kinds of edge information, DigraphEdge is a struct template.

template <typename EdgeInfo>
struct DigraphEdge
{
    int fromVertex;
    int toVertex;
    EdgeInfo einfo;
};



// A DigraphVertex includes two things: a VertexInfo object and a list of
// its outgoing edges.  Because different kinds of Digraphs store different
// kinds of vertex and edge information, DigraphVertex is a struct template.

template <typename VertexInfo, typename EdgeInfo>
struct DigraphVertex
{
    VertexInfo vinfo;
    std::list<DigraphEdge<EdgeInfo>> edges;
};



// Digraph is a class template that represents a directed graph implemented
// using adjacency lists.  It takes two type parameters:
//
// * VertexInfo, which specifies the kind of object stored for each vertex
// * EdgeInfo, which specifies the kind of object stored for each edge
//
// You'll need to implement the member functions declared here; each has a
// comment detailing how it is intended to work.
//
// Each vertex in a Digraph is identified uniquely by a "vertex number".
// Vertex numbers are not necessarily sequential and they are not necessarily
// zero- or one-based.

template <typename VertexInfo, typename EdgeInfo>
class Digraph
{
public:
    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    Digraph();

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    Digraph(const Digraph& d);

    // The move constructor initializes a new Digraph from an expiring one.
    Digraph(Digraph&& d) noexcept;

    // The destructor deallocates any memory associated with the Digraph.
    ~Digraph() noexcept;

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    Digraph& operator=(const Digraph& d);

    // The move assignment operator assigns the contents of an expiring
    // Digraph into "this" Digraph.
    Digraph& operator=(Digraph&& d) noexcept;

    // vertices() returns a std::vector containing the vertex numbers of
    // every vertex in this Digraph.
    std::vector<int> vertices() const;

    // edges() returns a std::vector of std::pairs, in which each pair
    // contains the "from" and "to" vertex numbers of an edge in this
    // Digraph.  All edges are included in the std::vector.
    std::vector<std::pair<int, int>> edges() const;

    // This overload of edges() returns a std::vector of std::pairs, in
    // which each pair contains the "from" and "to" vertex numbers of an
    // edge in this Digraph.  Only edges outgoing from the given vertex
    // number are included in the std::vector.  If the given vertex does
    // not exist, a DigraphException is thrown instead.
    std::vector<std::pair<int, int>> edges(int vertex) const;

    // vertexInfo() returns the VertexInfo object belonging to the vertex
    // with the given vertex number.  If that vertex does not exist, a
    // DigraphException is thrown instead.
    VertexInfo vertexInfo(int vertex) const;

    // edgeInfo() returns the EdgeInfo object belonging to the edge
    // with the given "from" and "to" vertex numbers.  If either of those
    // vertices does not exist *or* if the edge does not exist, a
    // DigraphException is thrown instead.
    EdgeInfo edgeInfo(int fromVertex, int toVertex) const;

    // addVertex() adds a vertex to the Digraph with the given vertex
    // number and VertexInfo object.  If there is already a vertex in
    // the graph with the given vertex number, a DigraphException is
    // thrown instead.
    void addVertex(int vertex, const VertexInfo& vinfo);

    // addEdge() adds an edge to the Digraph pointing from the given
    // "from" vertex number to the given "to" vertex number, and
    // associates with the given EdgeInfo object with it.  If one
    // of the vertices does not exist *or* if the same edge is already
    // present in the graph, a DigraphException is thrown instead.
    void addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo);

    // removeVertex() removes the vertex (and all of its incoming
    // and outgoing edges) with the given vertex number from the
    // Digraph.  If the vertex does not exist already, a DigraphException
    // is thrown instead. 
    void removeVertex(int vertex);

    // removeEdge() removes the edge pointing from the given "from"
    // vertex number to the given "to" vertex number from the Digraph.
    // If either of these vertices does not exist *or* if the edge
    // is not already present in the graph, a DigraphException is
    // thrown instead.
    void removeEdge(int fromVertex, int toVertex);

    // vertexCount() returns the number of vertices in the graph.
    int vertexCount() const noexcept;

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    int edgeCount() const noexcept;

    // This overload of edgeCount() returns the number of edges in
    // the graph that are outgoing from the given vertex number.
    // If the given vertex does not exist, a DigraphException is
    // thrown instead.
    int edgeCount(int vertex) const;

    // isStronglyConnected() returns true if the Digraph is strongly
    // connected (i.e., every vertex is reachable from every other),
    // false otherwise.
    bool isStronglyConnected() const;

    // findShortestPaths() takes a start vertex number and a function
    // that takes an EdgeInfo object and determines an edge weight.
    // It uses Dijkstra's Shortest Path Algorithm to determine the
    // shortest paths from the start vertex to every other vertex
    // in the graph.  The result is returned as a std::map<int, int>
    // where the keys are vertex numbers and the value associated
    // with each key k is the precedessor of that vertex chosen by
    // the algorithm.  For any vertex without a predecessor (e.g.,
    // a vertex that was never reached, or the start vertex itself),
    // the value is simply a copy of the key.
    std::map<int, int> findShortestPaths(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const;


private:
    // Add whatever member variables you think you need here.  One
    // possibility is a std::map where the keys are vertex numbers
    // and the values are DigraphVertex<VertexInfo, EdgeInfo> objects.


    // You can also feel free to add any additional member functions
    // you'd like (public or private), so long as you don't remove or
    // change the signatures of the ones that already exist.
    std::map<int, DigraphVertex<VertexInfo, EdgeInfo>> locations;

    bool edgeExist(int fromVertex, int toVertex) const;
    bool DFTr(int& start, std::map<int, bool>& visits, int& count) const;
    bool DFT(int& start) const;


};



// You'll need to implement the member functions below.  There's enough
// code in place to make them compile, but they'll all need to do the
// correct thing instead.

template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph()
{
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(const Digraph& d)
{
    // typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::iterator vertex; 
    // for(vertex = d.locations.begin(); vertex < d.locations.begin(); vertex++)
    // {
    //     locations.insert(std::pair(vertex->first, vertex->second));
    // }   

    locations = d.locations; 
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(Digraph&& d) noexcept
{
    std::swap(locations, d.locations);
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::~Digraph() noexcept
{
    locations.clear();
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(const Digraph& d)
{
    locations.clear();
    // typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::iterator vertex; 
    // for(vertex = d.locations.begin(); vertex < d.locations.begin(); vertex++)
    // {
    //     locations.insert(std::pair(vertex->first, vertex->second));
    // }   
    locations = d.locations; 

    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(Digraph&& d) noexcept
{
    locations.clear();
    std::swap(locations, d.locations);
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<int> Digraph<VertexInfo, EdgeInfo>::vertices() const
{   
    std::vector<int > vertexNumbers;
    for(auto it = locations.begin(); it != locations.end(); ++it )
    {
        vertexNumbers.push_back(it->first);
    }
    
    return vertexNumbers;

}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges() const
{
    std::vector<std::pair<int, int>> mapEdges;
    for(auto itL = locations.begin(); itL != locations.end(); ++itL)
    {
        for(auto it = locations.at(itL->first).edges.begin(); it != locations.at(itL->first).edges.end(); ++it)
        {   

            std::pair<int,int> temp;
            temp = std::make_pair(it->fromVertex, it->toVertex);
            mapEdges.push_back(temp);
        }

    }
    return mapEdges;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges(int vertex) const
{
    if(locations.count(vertex) == 0)
    {
        throw DigraphException("Vertex does not exist");
    }

    //throw DigraphException("Vertex does not exist");
    std::vector<std::pair<int, int>> vertexEdges;
    for(auto it = locations.at(vertex).edges.begin(); it!= locations.at(vertex).edges.end(); it++)
    {
        std::pair<int,int> temp;
        temp = std::make_pair(it->fromVertex, it->toVertex);
        vertexEdges.push_back(temp);        
    }
    return vertexEdges;
}


template <typename VertexInfo, typename EdgeInfo>
VertexInfo Digraph<VertexInfo, EdgeInfo>::vertexInfo(int vertex) const
{
    if(locations.find(vertex) == locations.end())
    {
        throw DigraphException("Vertex does not exist");
    }
    VertexInfo temp;
    temp = locations.at(vertex).vinfo;
    return temp;
}


template <typename VertexInfo, typename EdgeInfo>
EdgeInfo Digraph<VertexInfo, EdgeInfo>::edgeInfo(int fromVertex, int toVertex) const
{

    EdgeInfo temp;
    if(locations.find(fromVertex) == locations.end())
    {
        throw DigraphException("Starting Vertex does not exist");
    }
    if(locations.find(toVertex) == locations.end())
    {
        throw DigraphException("End vertex does not exist");
    }
    if(edgeExist(fromVertex, toVertex)==false)
    {
        throw DigraphException("Edge does not exist");
    }

    for(auto it = locations.at(fromVertex).edges.begin(); it != locations.at(fromVertex).edges.end(); ++it)
    {
        if(it->toVertex == toVertex)
        {
            temp = it->einfo;
        }
    }
    
    return temp;

    // for(auto& it:locations.at(fromVertex).edges)
    // {
    //     if(it->toVertex == toVertex)
    //     {
    //         EdgeInfo temp;
    //         temp = it->einfo;
    //         return temp;
    //     }
    // }


}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addVertex(int vertex, const VertexInfo& vinfo)
{
    if(locations.size() == 0)
    {
        DigraphVertex<VertexInfo, EdgeInfo> current;
        current.vinfo = vinfo;
        locations[vertex] = current;

    }
    else
    {
        for(auto& it:locations)
        {
            if(vertex == it.first)
            {
                throw DigraphException("Vertex already exists");
            }
        }

        DigraphVertex<VertexInfo, EdgeInfo> current;
        current.vinfo = vinfo;
        locations[vertex] = current;
    }
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo)
{

    if(locations.find(fromVertex) == locations.end())
    {
        throw DigraphException("Starting Vertex does not exist");
    }
    if(locations.find(toVertex) == locations.end())
    {
        throw DigraphException("End vertex does not exist");
    }

    if(locations[fromVertex].edges.size() == 0)
    {
        DigraphEdge<EdgeInfo> temp;
        temp.fromVertex = fromVertex;
        temp.toVertex = toVertex;
        temp.einfo = einfo;
        auto it = locations[fromVertex].edges.begin();
        locations[fromVertex].edges.insert(it, temp);
    }
    else
    {
        for(auto it = locations[fromVertex].edges.begin(); it != locations[fromVertex].edges.end(); ++it)
        {
            if(it->toVertex == toVertex)
            {
                throw DigraphException("Edge already exists");
            }
        }
        DigraphEdge<EdgeInfo> temp;
        temp.fromVertex = fromVertex;
        temp.toVertex = toVertex;
        temp.einfo = einfo;
        auto it = locations[fromVertex].edges.begin();
        locations[fromVertex].edges.insert(it, temp);
    }   
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeVertex(int vertex)
{
    if(locations.find(vertex) == locations.end())
    {
        throw DigraphException("Vertex does not exist");
    }
    
    locations.erase(vertex);
    for(auto it = locations.begin(); it != locations.end();++it)
    {
        for(auto it2 = locations.at(it->first).edges.begin(); it2 != locations.at(it->first).edges.end(); ++it2)
        {
            if(it2 -> toVertex == vertex)
            {
                locations.at(it->first).edges.erase(it2);
            }
        }
    }

    
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeEdge(int fromVertex, int toVertex)
{
    if(locations.find(fromVertex) == locations.end())
    {
        throw DigraphException("Starting Vertex does not exist");
    }
    if(locations.find(toVertex) == locations.end())
    {
        throw DigraphException("End vertex does not exist");
    }
    if(edgeExist(fromVertex, toVertex) == false)
    {
        throw DigraphException("Edge does not exist");
    }

    for(auto it = locations.at(fromVertex).edges.begin(); it != locations.at(fromVertex).edges.end(); ++it)
    {
        if(it->toVertex == toVertex)
        {
            locations.at(fromVertex).edges.erase(it);
        }
    }

}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::vertexCount() const noexcept
{
    return locations.size();
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount() const noexcept
{
    int numEdges = 0;
    for(auto it = locations.begin(); it != locations.end(); ++it)
    {
        numEdges += locations.at(it->first).edges.size();
    }
    return numEdges;


}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount(int vertex) const
{
    
    if(locations.count(vertex)==0)
    {
        throw DigraphException("Vertex does not exist");
    }
    int numEdges = locations.at(vertex).edges.size();
    
    return numEdges;
}


template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::isStronglyConnected() const
{
    bool connected = false;
    std::vector<int> vertices = this->vertices();
    if(vertices.size() == 0)
    {
        return true;
    }
    int numVertices = vertices.size();

    for(int i =0; i<numVertices; ++i)
    {
       connected = DFT(vertices.at(i));
        //std::cout << "calling DFT on " << vertices.at(i) << "connected is " << connected <<std::endl;
    }

    return connected;
}


template <typename VertexInfo, typename EdgeInfo>
std::map<int, int> Digraph<VertexInfo, EdgeInfo>::findShortestPaths(
    int startVertex,
    std::function<double(const EdgeInfo&)> edgeWeightFunc) const
{
    // map for final paths. Pv predecessor
    std::map<int, int> shortestPaths;
    //map for Kv(shortest path known)and Dv(shortest distance found so far)
    std::map<int, std::pair<bool, double>> info;
    for(auto it = locations.begin(); it!=locations.end();++it)
    {
        shortestPaths.insert(std::make_pair(it->first, -99));
        info.insert(std::make_pair(it->first, std::make_pair(false, 9999999)));
    }

    shortestPaths[startVertex] = startVertex;
    info[startVertex].second = 0;

    //min heap of pairs<weight, vertex>;
    typedef std::pair<int, int> P;
    std::priority_queue< P, std::vector<P>, std::greater<P> > pq;
    pq.push(std::make_pair(0, startVertex));
    while(pq.empty() == false)
    {
        int vertex = pq.top().second;
        pq.pop();
        if(info[vertex].first == false)
        {
            info[vertex].first = true;
            for(auto it = locations.at(vertex).edges.begin(); it != locations.at(vertex).edges.end(); ++it)
            {
                //new path length
                double distance = info[vertex].second + edgeWeightFunc(it->einfo); 
                if(info[it->toVertex].second > distance)
                {
                    //std::cout << info[it->toVertex].second << std::endl;
                    info[it->toVertex].second = distance;
                    shortestPaths[it->toVertex] = vertex;
                    //std::cout << it->toVertex << "'s Pred is Vertex, and their distance is " << info[it->toVertex].second << std::endl;
                    pq.push(std::make_pair(distance, it->toVertex));
                }
                //std::cout << it->toVertex << " " << info[vertex].second << " " << distance << std::endl;
            }
        }
    }

    return shortestPaths;
}

template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::edgeExist(int fromVertex, int toVertex) const
{
    bool exist = false; 
    for(auto it = locations.at(fromVertex).edges.begin(); it != locations.at(fromVertex).edges.end();++it)
    {
        if(it->toVertex == toVertex)
        {
            exist = true;
        }
    }
    return exist;
}

template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::DFTr(int& start, std::map<int, bool>& visits, int& count) const
{
    bool connected = false;
    std::vector<std::pair<int, int>> temp = this->edges(start);
    for(int i=0; i<temp.size();++i)
    {
        if(visits.at(temp.at(i).second) == false)
        {
            visits.at(temp.at(i).second) = true;
            count++;
            //std::cout << "Visited " << temp.at(i).second << " Count is now " << count << std::endl;
            if(count == visits.size())
            {
                
                connected = true;
               // std::cout << "vertex " <<  temp.at(i).second << " links to all other vertices, connected is "<< connected << std::endl;
                return connected;
            }
            connected = DFTr(temp.at(i).second, visits, count);
        }
    }

    // if(count == visits.size())
    // {           
    //     connected = true;
    //     //std::cout << "vertex " <<  temp.at(i).second << " links to all other vertices, connected is "<< connected << std::endl;
    //     return connected;
    // }    
    //std::cout << " some how it got to here and connecte is now" << connected << "and count is now " << count << std::endl;
    return connected;

}

template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::DFT(int& start) const
{
    bool connected = false;
    int visitCount = 0;
    std::vector<int> vertices = this->vertices();
    int numVertices = vertices.size();
    std::map<int, bool> visits;
    for(int i =0; i<numVertices; ++i)
    { 
        visits.insert(std::pair(vertices.at(i), false));
    }

    connected = DFTr(start, visits, visitCount);
    //std::cout << " connected is " << connected << "after returning to DFT" << std::endl;

    return connected; 
}




#endif // DIGRAPH_HPP

