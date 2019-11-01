#include "network.h"
#include "random.h"
#include <map>

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}


std::pair<size_t, double> Network::degree(const size_t& n) const{
	double sum_intens;
	size_t nbr_connec(neighbors(n).size());
	
	for( auto i : neighbors(n)){
		sum_intens += i.second;
		}
		
	return std::make_pair(nbr_connec, sum_intens);
}

//trouve les voisons du neuron en paramètre
std::vector<std::pair<size_t, double> > Network::neighbors(const size_t& nrn) const{
	std::vector<std::pair<size_t, double> > vec;
	std::map<std::pair<size_t, size_t>,double >::const_iterator it;	//const_iterator car la méthode est const donc on peut pas autoriser qqch qui peut modifier

	for ( it = links.lower_bound({nrn,0}); it !=links.cend() and nrn == it->first.first; ++it){
		vec.push_back(std::make_pair(it->first.second, it->second));
	}
	return vec;
}


std::set<size_t> Network::step(const std::vector<double>& vec){
	std::set<size_t> set;
	for (size_t m(0); m<neurons.size(); ++m){ 
		if (neurons[m].firing()) { 
			neurons[m].reset();	//il faut reset les firing neurons avant tout.
			set.insert(m);
		}		
	}

	for ( size_t i(0) ; i< neurons.size(); ++i) {
		double excit_sum(0.0);
		double inhib_sum(0.0);
		
		double w(1);
		if (neurons[i].is_inhibitory()) { w = 0.4;}
		
		
		std::vector<std::pair<size_t, double>> sending_neurons = neighbors(i); 
		
		for (size_t j(0) ; j< sending_neurons.size(); ++j) {
			size_t n = sending_neurons[j].first;	//trouve le jème voisin du neuron i et donne son indice à n
			
			if (set.count(n) != 0){
				if(neurons[n].is_inhibitory()) {
					inhib_sum += neighbors(i)[j].second;
				}
				else { 
					excit_sum += neighbors(i)[j].second;
				}
			} 
		}
		
		double intensity = w * vec[i] + 0.5*excit_sum + inhib_sum; 
		neurons[i].input(intensity);      
		neurons[i].step();
	}
	
	return set;
}


