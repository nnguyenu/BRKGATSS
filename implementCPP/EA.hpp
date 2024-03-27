#include "utils.hpp"

std::vector<int> randomBitstring(int n){
    srand(time(NULL));
    std::vector<int> bitstring(n);
    for(int i = 0; i < n; i++){
        bitstring[i] = rand() % 2;
    }
    return bitstring;
}

std::vector<int> randomBitstringR(int n, double rate){
    srand(time(NULL));
    std::vector<int> bitstring(n, 0);
    for(int i = 0; i < n; i++){
        if(randomFloat() < rate){
            bitstring[i] = 1;
        }
    }
    return bitstring;
}

class EA {
    public:
        std::vector<std::vector<int>> G;
        int n, m;
        double time_limit = 100;
        std::vector<int> deg;
        std::vector<int> bitstring;
        std::vector<std::pair<int,int>> deg_id;
        EA() = default;
        EA(std::vector<std::vector<int>>graph, int v,int e){
            G = graph;
            n = v;
            m = e;
            deg.resize(n);
            deg_id.resize(n);
            for(int i = 0; i < n; i++){
                deg[i] = G[i].size();
                deg_id[i] = {deg[i], i};
            }
            srand(time(NULL));
            std::sort(deg_id.begin(), deg_id.end());
            bitstring = std::vector<int>(n, 0);
            time_limit = std::max((double)100, (double)n / 100);
        }
        double get_mutation_rate(){
            return 1.0 / n;
        }
        std::vector<int> mutate(){
            double mutation_rate = get_mutation_rate();
            std::vector<int> newbitstring = bitstring;
            for(int i=0;i<n;i++){
                double x = randomFloat();
                if(x < mutation_rate){
                    newbitstring[i] = 1 - newbitstring[i];
                }
            }
            return newbitstring;
        }
        std::vector<int> MDG(std::vector<int>S){
            auto [Cov, cnt] = phiBig(S);
            int sum = std::accumulate(Cov.begin(), Cov.end(), 0);
            if(sum == n)    return S;
            for(int i=deg_id.size()-1;i>=0;i--){
                int id = deg_id[i].second;
                if(Cov[id] == 0){
                    S[id] = 1;
                    std::tie(Cov, cnt) = phi(Cov, cnt, id);
                    sum = std::accumulate(Cov.begin(), Cov.end(), 0);
                    if(sum == n)    break;
                }
            }
            return S;
        }
        
        // Complexity: O(m + n)
        std::pair<std::vector<int>,std::vector<int>> phiBig(std::vector<int> Cov){
            std::queue<int> q;
            for(int i = 0; i < n; i++) if(Cov[i] == 1) q.push(i);  
            std::vector<int> d(n, 0), res(n, 0);    
            while(!q.empty()){  
                int u = q.front();q.pop();
                if(res[u] == 1)     continue;
                res[u] = 1;
                for(auto v : G[u]){
                    d[v] += 1;
                    if(d[v] == (deg[v]+1)/2)
                        q.push(v);
                }
            }
            return {res, d};
        }
        std::pair<std::vector<int>,std::vector<int>> phi(std::vector<int> Cov, std::vector<int> d, int v){
            std::queue<int> q;q.push(v); 
            std::vector<int> res = Cov;
            while(!q.empty()){  
                int u = q.front();q.pop();
                if(res[u] == 1)     continue;
                res[u] = 1;
                for(auto v : G[u]){
                    d[v] += 1;
                    if(d[v] == (deg[v]+1)/2)
                        q.push(v);
                }
            }
            return {res, d};
        }
        // Complexity: n * (m + n)
        std::vector<int> reverseMDG(std::vector<int> S){
            //std::cout << "Previous: " << std::accumulate(S.begin(), S.end(), 0);
            for(auto [_, i] : deg_id){
                if(S[i] == 1){
                    S[i] = 0;
                    auto [eval, _] = phiBig(S);
                    if(std::accumulate(eval.begin(), eval.end(), 0) != n){
                        S[i] = 1;
                    }
                }
            }
            //std::cout << " After: " << std::accumulate(S.begin(), S.end(), 0) << '\n';
            return S;
        }

        std::pair<int, std::vector<int>> fitness(std::vector<int> S){
            std::vector<int> mdg = MDG(S);
            mdg = reverseMDG(mdg);
            int sum = std::accumulate(mdg.begin(), mdg.end(), 0);
            return {n - sum, mdg};
        }

        void run(int threshold = 0, std::string filename = "EA.txt"){
            // let x be the value and vector xRes be the result of fitness function
            int x, y;
            std::vector<int> xRes, yRes;
            std::tie(x, xRes) = fitness(bitstring);
            //std::cout << "Initial: " << n - x << '\n';
            auto start = std::chrono::high_resolution_clock::now();
            double besttime;
            int cnt = 0, bestcnt;
            while(1){
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                if(duration.count() > time_limit)   break;
                std::vector<int> newbitstring = mutate();
                std::tie(y, yRes) = fitness(newbitstring);
                
                if(x < y){
                    bitstring = newbitstring;
                    x = y;
                    xRes = yRes;
                    besttime = duration.count();
                    bestcnt = cnt;
                    
                    std::ofstream outfile;outfile.open(filename, std::ios_base::app);
                    outfile << duration.count() << ' ' << cnt << ' ' << n - x << ' ';
                    for(int i=0;i<n;i++)    if(xRes[i] == 1)   outfile << i << ' ';
                    outfile << '\n';
                    outfile.close();
                    
                    //std::cout << ">" << duration.count() << " New best: " << n - x << '\n';
                    if(n - x <= threshold)  break;
                }
                cnt++;
            }
            std::ofstream outfile;outfile.open(filename, std::ios_base::app);
            
            outfile << "Finish " << besttime << ' ' << n - x << ' ' << bestcnt << ' ';
            for(int i=0;i<n;i++)    if(xRes[i] == 1)   outfile << i << ' ';
            outfile << '\n';
            
            outfile.close();
            // std::cout << "Best: " << n - x << " Time: " << besttime << '\n';
            // std::cout << cnt << '\n';
        }
};

class fastEA : public EA {
    public:
        double beta = 1.5;
        std::vector<double> pm_distribution, pm_power_law;
        fastEA(std::vector<std::vector<int>>graph, int v,int e): EA(graph, v, e){
            std::vector<double> prefix_sum;
            pm_distribution.resize(n/2);
            prefix_sum.resize(pm_distribution.size());
            for(int i=0;i<pm_distribution.size();i++){
                pm_distribution[i] = pow(i+1, -beta);
                prefix_sum[i] = pm_distribution[i];
                if(i != 0)  prefix_sum[i] += prefix_sum[i-1];
            }    
        }
        double get_mutation_rate(){
            int id = sample(pm_distribution);
            return id / n;
        }
};

class balancedEA : public EA {
    public:
        std::random_device rd;     // Only used once to initialise (seed) engine
        std::mt19937 rng;    
        std::uniform_int_distribution <int> uni;
        balancedEA(std::vector<std::vector<int>>graph, int v,int e): EA(graph, v, e){
            rng = std::mt19937(rd());    // Random-number engine used (Mersenne-Twister in this case)
            uni = std::uniform_int_distribution<int>(0, n-1); // guaranteed unbiased
        }
        std::vector<int> mutate(){
            double p = randomFloat();
            std::vector<int> newbitstring = bitstring;
            if(p > 0.5){
                double mutation_rate = get_mutation_rate();
                for(int i=0;i<n;i++){
                    double x = randomFloat();
                    if(x < mutation_rate){
                        newbitstring[i] = 1 - newbitstring[i];
                    }
                }
                return newbitstring;
            }else{
                // choose random vertex from 0 to n-1
                int v = uni(rng);
                std::vector<int> Nv;
                for(auto i: G[v]){
                    if(bitstring[i] != bitstring[v]){
                        Nv.push_back(i);
                    }
                }
                if(Nv.empty())  return bitstring;
                int u = Nv[ rand() % Nv.size()];
                newbitstring[v] = 1 - newbitstring[v];
                newbitstring[u] = 1 - newbitstring[u];
                return newbitstring;
            }
        }
};