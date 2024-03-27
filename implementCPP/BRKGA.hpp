#include "utils.hpp"

class BRKGA {
    public:
        std::vector<std::vector<int>> G;
        int n, m;
        double time_limit = 100;
        int n_ind = 46;
        double pe = 0.24, pm = 0.13, pelite = 0.69;
        bool seed = true;
        std::vector<int> deg;
        std::vector<int> score;
        std::vector<std::vector<double>> P;
        std::vector<std::vector<double>> Pm;
        std::vector<std::vector<double>> Pc;
        std::vector<std::vector<int>> result;
        int Pm_size, Pe_size, Pc_size;
        
        BRKGA() = default;
        BRKGA(std::vector<std::vector<int>>graph, int v,int e){
            G = graph;
            n = v;
            m = e;
            deg.resize(n);
            for(int i = 0; i < n; i++){
                deg[i] = G[i].size();
            }
            score.resize(n_ind);
            P.resize(n_ind);
            Pm.resize(n_ind);
            Pc.resize(n_ind);
            result.resize(n_ind);
            for(int i=0;i<n_ind;i++){
                P[i].resize(n);
                Pm[i].resize(n);
                Pc[i].resize(n);
                result[i].resize(n);
            }
            time_limit = std::max((double)100, (double)n / 100);
            srand(time(0));
            for(int i = 0; i < n_ind; i++){
                for(int j = 0; j < n; j++){
                    P[i][j] = randomFloat();
                }
            }
            if(seed == true)    P[0] = std::vector<double>(n, 0.5);
        }

        double get_pe(){
            return pe;
        }
        double get_pm(){
            return pm;
        }
        double get_pelite(){
            return pelite;
        }
        std::vector<int> run(int threshold = 0, std::string filename="BRKGA.txt"){
            srand(time(0));
            for(int i = 0; i < n_ind; i++){
                for(int j = 0; j < n; j++){
                    P[i][j] = randomFloat();
                }
            }
            int cnt = 0;
            int best = n, bestcnt;
            if(seed == true)    P[0] = std::vector<double>(n, 0.5);
            auto start = std::chrono::high_resolution_clock::now();
            double besttime;
            std::vector<int> bestres;
            while(1){
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                if(duration.count() > time_limit)   break;

                std::pair<std::vector<int>, std::vector<int>> id = get_elite_and_nonelite_id();
                std::vector<int> elite_id = id.first;
                std::vector<int> non_elite_id = id.second;
                mutate();
                Pc_size = n_ind - Pm_size - elite_id.size();
                crossover(elite_id, non_elite_id);
                std::vector<std::vector<double>> newPn(elite_id.size());
                for(int i = 0; i < elite_id.size(); i++)    {newPn[i] = P[elite_id[i]];}
                newPn.insert(newPn.end(), Pm.begin(), Pm.begin() + Pm_size);
                newPn.insert(newPn.end(), Pc.begin(), Pc.begin() + Pc_size);
                P = newPn;

                int bestVal = *std::min_element(score.begin(), score.end());
                if(best > bestVal){
                    best = bestVal;
                    besttime = duration.count();
                    bestcnt = cnt;
                    int bestValId = std::find(score.begin(), score.end(), bestVal) - score.begin();
                    bestres = result[bestValId];
                    std::ofstream outfile;outfile.open(filename, std::ios_base::app);
                    outfile << duration.count() << ' ' << cnt << ' ' << best << ' ';
                    for(int i=0;i<n;i++)    if(result[bestValId][i] == 1)   outfile << i << ' ';
                    outfile << '\n';
                    outfile.close();
                    //std::cout << '>' << duration.count() << " New best: " << best << '\n';
                    if(best <= threshold)   break;
                }
                cnt++;
            }
            std::ofstream outfile;outfile.open(filename, std::ios_base::app);
            outfile << "Finish " << besttime << ' ' << best << ' ' << bestcnt << ' ';
            for(int i=0;i<n;i++)    if(bestres[i] == 1)   outfile << i << ' ';
            outfile << '\n';
            outfile.close();
            std::cout << "Best: " << best << " Time: " << besttime << '\n';

            return bestres;
        }

        std::pair<std::vector<int>, std::vector<int>> get_elite_and_nonelite_id(){
            double cur_pelite = get_pe();
            Pe_size = ceil(cur_pelite * n_ind);
            eval();
            std::vector<std::pair<int,int>> eid(n_ind);
            for(int i = 0; i < n_ind; i++){
                eid[i] = std::make_pair(score[i], i);
            }
            std::sort(eid.begin(), eid.end());
            // return the 2nd value for the first Pe_size elements
            std::vector<int> res1(Pe_size), res2(n_ind - Pe_size);
            for(int i = 0; i < Pe_size; i++){
                res1[i] = eid[i].second;
            }
            for(int i = Pe_size; i < n_ind; i++){
                res2[i - Pe_size] = eid[i].second;
            }
            return std::make_pair(res1, res2);
        }

        void mutate(){
            double cur_pm = get_pm();
            Pm_size = ceil(cur_pm * n_ind);
            for(int i = 0; i < Pm_size; i++){
                for(int j = 0; j < n; j++){
                    Pm[i][j] = randomFloat();
                }
            }
        }

        void crossover(std::vector<int> elite, std::vector<int> non_elite){
            double cur_pe = get_pelite();
            for(int i=0;i<Pc_size;i++){
                int x = *select_randomly(elite.begin(), elite.end());
                int y = *select_randomly(non_elite.begin(), non_elite.end());
                for(int j=0;j<n;j++){
                    if(randomFloat() < cur_pe)  Pc[i][j] = P[x][j];
                    else                        Pc[i][j] = P[y][j];
                    
                }
            }
        }
       
        void eval(){
            std::vector<double> d(n);
            for(int i = 0; i < n_ind; i++){
                for(int j=0;j<n;j++){
                    d[j] = deg[j] * P[i][j];
                }
                result[i] = MDG(d);
                int x = std::accumulate(result[i].begin(), result[i].end(), 0);
                score[i] = x;
            }
        }
        // Complexity: MDG * (n+m) with MDG be the solution value
        std::vector<int> MDG(std::vector<double> d){
            std::vector<int> Cov(n, 0), S(n, 0);
            while(1){
                int sum = std::accumulate(Cov.begin(), Cov.end(), 0);
                if(sum == n)    break;
                int maxd = -1, id = -1;
                for(int i=0;i<n;i++){
                    if(Cov[i] == 0){
                        if(d[i] > maxd){
                            maxd = d[i];
                            id = i;
                        }
                    }
                }
                Cov[id] = 1;
                S[id] = 1;
                Cov = phi(Cov);
            }
            return S;
        }
        // Complexity: O(m + n)
        std::vector<int> phi(std::vector<int> Cov){
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
            return res;
        }
};

std::ostream& operator<< (std::ostream &out, BRKGA const& data) {
    out << "n: " << data.n << " m: " << data.m;
    return out;
}

class fastBRKGA: public BRKGA {
    public:
        double beta = 1.5;
        std::vector<double> pe_distribution, pe_power_law;
        std::vector<double> pm_distribution, pm_power_law;
        std::vector<double> pelite_distribution, pelite_power_law;

        // inherit everything from BRKGA
        fastBRKGA(std::vector<std::vector<int>>graph, int v,int e) : BRKGA(graph, v, e)
        {
            std::vector<double> prefix_sum;

            pe_distribution.resize(16);
            prefix_sum.resize(pe_distribution.size());
            for(int i=0;i<pe_distribution.size();i++){
                pe_distribution[i] = pow(i+1, -beta);
                prefix_sum[i] = pe_distribution[i];
                if(i != 0)  prefix_sum[i] += prefix_sum[i-1];
            }
            pe_power_law.resize(16);
            for(int i=0;i<pe_power_law.size();i++){
                pe_power_law[i] = prefix_sum[i] / prefix_sum.back();
            }

            pm_distribution.resize(20);
            prefix_sum.resize(pm_distribution.size());
            for(int i=0;i<pm_distribution.size();i++){
                pm_distribution[i] = pow(i+1, -beta);
                prefix_sum[i] = pm_distribution[i];
                if(i != 0)  prefix_sum[i] += prefix_sum[i-1];
            }
            pm_power_law.resize(20);
            for(int i=0;i<pm_power_law.size();i++){
                pm_power_law[i] = prefix_sum[i] / prefix_sum.back();
            }

            pelite_distribution.resize(30);
            prefix_sum.resize(pelite_distribution.size());
            for(int i=0;i<pelite_distribution.size();i++){
                pelite_distribution[i] = pow(i+1, -beta);
                prefix_sum[i] = pelite_distribution[i];
                if(i != 0)  prefix_sum[i] += prefix_sum[i-1];
            }
            pelite_power_law.resize(30);
            for(int i=0;i<pelite_power_law.size();i++){
                pelite_power_law[i] = prefix_sum[i] / prefix_sum.back();
            }
        }
        double get_pe(){
            int x = sample(pe_power_law);
            return 0.1 + 0.01 * (15 - x);
        }
        double get_pm(){
            int x = sample(pm_power_law);
            return 0.1 + 0.01 * x;
        }
        double get_pelite(){
            int x = sample(pelite_power_law);
            return 0.5 + 0.01 * x;
        }
};