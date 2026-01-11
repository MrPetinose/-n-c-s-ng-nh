/**
 * A Trainable Clustering Algorithm based on Shortest Paths from Density Peaks
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <queue>
#include <cstdlib>
#include <chrono>

using namespace std;

// Cấu trúc dữ liệu cho một điểm
struct Point {
    int id;
    double x, y;
    double rho;         // Mật độ
    double delta;       // Khoảng cách đến điểm đặc hơn
    double gamma;       // Gamma = Rho * Delta
    int nearestHigher;  // ID của điểm đặc hơn gần nhất
    int clusterId;      // Nhãn cụm
    bool isPeak;
};

// Cấu trúc dùng cho Dijkstra Priority Queue
struct NodeState {
    int u;
    double cost;
    bool operator>(const NodeState& other) const {
        return cost > other.cost;
    }
};

class AutoClustering {
private:
    vector<Point> points;
    vector<vector<double>> distMatrix;
    int N;
    double dc; // Khoảng cách cắt tự động

public:
    // 1. Đọc dữ liệu
    void loadData(string filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "[Loi] Khong the mo file " << filename << endl;
            exit(1);
        }

        string line;
        // Bỏ qua dòng tiêu đề nếu có (chứa chữ cái)
        getline(file, line);
        if (line.find_first_not_of("0123456789. \t") == string::npos) {
            // Nếu dòng đầu toàn số, reset lại để đọc từ đầu
            file.clear();
            file.seekg(0);
        }

        while (getline(file, line)) {
            if (line.empty()) continue;
            stringstream ss(line);
            Point p;
            ss >> p.id >> p.x >> p.y;
            
            // Reset các chỉ số
            p.rho = 0; p.delta = 0; p.gamma = 0;
            p.nearestHigher = -1; p.clusterId = -1; p.isPeak = false;
            
            points.push_back(p);
        }
        N = points.size();
        if (N < 2) {
            cerr << "[Loi] Du lieu qua it (N < 2). Khong the phan cum." << endl;
            exit(1);
        }
        cout << "-> Da tai " << N << " diem du lieu." << endl;
        file.close();
    }

    // 2. Tính Ma trận khoảng cách & Tự động chọn dc 
    void calculateDistanceMatrix() {
        distMatrix.resize(N, vector<double>(N));
        vector<double> sampleDistances;
        sampleDistances.reserve(N * N / 2);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    distMatrix[i][j] = 0;
                } else {
                    double dx = points[i].x - points[j].x;
                    double dy = points[i].y - points[j].y;
                    distMatrix[i][j] = sqrt(dx*dx + dy*dy);
                    
                    if (i < j) {
                        sampleDistances.push_back(distMatrix[i][j]);
                    }
                }
            }
        }

        // Sắp xếp để tìm Percentile
    sort(sampleDistances.begin(), sampleDistances.end());
        int idx;
        if (sampleDistances.size() < 1000) {
            idx = (int)(sampleDistances.size() * 0.20); 
        } else {
            idx = (int)(sampleDistances.size() * 0.02);
        }
        idx = max(0, idx); 
        dc = sampleDistances[idx];

        if (dc <= 1e-9) { 
            int fallbackIdx = max(0, (int)(sampleDistances.size() * 0.10));
            dc = sampleDistances[fallbackIdx];
        }
        if (dc <= 1e-9) dc = 1.0; 

        cout << "-> Cutoff Distance: " << dc << endl;
    }

    // 3. Tính Mật độ
    void calculateRho() {
        for (int i = 0; i < N; ++i) {
            points[i].rho = 0;
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                points[i].rho += exp(-pow(distMatrix[i][j] / dc, 2));
            }
        }
    }

    // 4. Tính Delta
    void calculateDelta() {
        for (int i = 0; i < N; ++i) {
            points[i].delta = numeric_limits<double>::max();
            points[i].nearestHigher = -1;
            
            bool foundHigher = false;
            for (int j = 0; j < N; ++j) {
                if (points[j].rho > points[i].rho) {
                    if (distMatrix[i][j] < points[i].delta) {
                        points[i].delta = distMatrix[i][j];
                        points[i].nearestHigher = j;
                        foundHigher = true;
                    }
                }
            }
            // Xử lý ngoại lệ cho điểm có mật độ cao nhất toàn cục
            if (!foundHigher) {
                double maxDist = 0;
                for(int j=0; j<N; ++j) {
                    if (distMatrix[i][j] > maxDist) maxDist = distMatrix[i][j];
                }
                points[i].delta = maxDist;
            }
            
            // Gamma của điểm 
            points[i].gamma = points[i].rho * points[i].delta;
        }
    }

    // 5. TÌM ĐỈNH
    vector<int> findPeaksAuto() {
        vector<int> peaks;
        
        // Tính trung bình và độ lệch chuẩn của Gamma
        double sumGamma = 0, sumSqGamma = 0;
        for (const auto& p : points) {
            sumGamma += p.gamma;
            sumSqGamma += p.gamma * p.gamma;
        }
        double mean = sumGamma / N;
        double variance = (sumSqGamma / N) - (mean * mean);
        double stdDev = sqrt(variance);

        cout << "-> Thong ke Gamma: Mean=" << mean << ", StdDev=" << stdDev << endl;

        double threshold = mean + 2.0 * stdDev; 

        if (stdDev < 1e-5) threshold = mean; 

        cout << "-> Nguong tu dong (Threshold): " << threshold << endl;
        cout << "-> Cac Dinh duoc chon:" << endl;
        cout << setw(5) << "ID" << setw(10) << "Gamma" << setw(10) << "Rho" << setw(10) << "Delta" << endl;

        // Chọn đỉnh
        vector<pair<double, int>> candidates;
        for(int i=0; i<N; ++i) candidates.push_back({points[i].gamma, i});
        sort(candidates.rbegin(), candidates.rend());

        for (int i = 0; i < N; ++i) {
            int idx = candidates[i].second;
            double val = candidates[i].first;

            // Điều kiện chọn: Gamma > Threshold
            if (val > threshold || peaks.empty()) {
                peaks.push_back(idx);
                points[idx].isPeak = true;
                points[idx].clusterId = points[idx].id;
                cout << setw(5) << points[idx].id 
                     << setw(10) << fixed << setprecision(2) << val
                     << setw(10) << points[idx].rho 
                     << setw(10) << points[idx].delta << endl;
            } else {
                if (peaks.size() >= 1) break; 
            }
        }     
        // if (peaks.size() > N / 5 && N > 50) {
        //     cout << "[Info] Phat hien qua nhieu dinh, giu lai Top 5..." << endl;
        //     peaks.resize(5);
        // }
        return peaks;
    }

    // 6. Hàm chi phí cạnh (Minimax - Euclidean)
    double getEdgeCost(int u, int v) {
        // [TRAINABLE PLACEHOLDER]
        return distMatrix[u][v];
    }

    // 7. Chạy Dijkstra Minimax (Tối ưu hóa toàn cục)
    void runGlobalOptimization(vector<int>& peaks) {
        vector<double> minMaxCost(N, numeric_limits<double>::max());
        vector<int> rootPeak(N, -1);
        priority_queue<NodeState, vector<NodeState>, greater<NodeState>> pq;

        // Khởi tạo nguồn ảo
        for (int peakIdx : peaks) {
            minMaxCost[peakIdx] = 0;
            rootPeak[peakIdx] = peakIdx;
            pq.push({peakIdx, 0});
        }

        while (!pq.empty()) {
            NodeState current = pq.top();
            pq.pop();

            int u = current.u;
            if (current.cost > minMaxCost[u]) continue;

            for (int v = 0; v < N; ++v) {
                if (u == v) continue;

                double edgeWeight = getEdgeCost(u, v);
                double newCost = max(current.cost, edgeWeight); // Minimax Logic

                if (newCost < minMaxCost[v]) {
                    minMaxCost[v] = newCost;
                    rootPeak[v] = rootPeak[u];
                    pq.push({v, newCost});
                }
            }
        }

        // Gán nhãn
        for (int i = 0; i < N; ++i) {
            if (rootPeak[i] != -1)
                points[i].clusterId = points[rootPeak[i]].id;
        }
    }

    // Xuất kết quả
    void exportResults(string filename) {
        ofstream file(filename);
        file << "Id,X,Y,Cluster_ID,Rho,Delta,Gamma,Is_Peak" << endl;
        for (const auto& p : points) {
            file << p.id << "," << fixed << setprecision(4) 
                 << p.x << "," << p.y << "," 
                 << p.clusterId << "," 
                 << setprecision(2) << p.rho << "," 
                 << p.delta << ","
                 << p.gamma << ","
                 << (p.isPeak ? "YES" : "NO") << endl;
        }
        file.close();
    }
};

int main() {
    AutoClustering ac;

    string inputFile = "../data/Random200Points.txt";
    //string inputFile = "../data/random_1000_points.txt";
    //string inputFile = "../data/data_2000_points.txt";
    //string inputFile = "../data/data_2000_complex_with_noise.txt";

    auto start = std::chrono::high_resolution_clock::now();

    cout << "--- BAT DAU  ---" << endl;
    ac.loadData(inputFile);
    ac.calculateDistanceMatrix();
    ac.calculateRho();
    ac.calculateDelta();
    vector<int> peaks = ac.findPeaksAuto();
    ac.runGlobalOptimization(peaks);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Tong thoi gian chay thuat toan: " << duration.count() << " ms (mili-giay)" << endl;

    cout << "--- HOAN THANH ---" << endl;
    
    ac.exportResults("Auto_Result.csv");
    int result = system("python Visualize_Result2.py");
    return 0;
}