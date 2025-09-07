#include "solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_map>
#include <vector>

using namespace std;

static vector<vector<double>> dist_btw_villages;
static vector<vector<double>> dist_btw_cities_villages;

struct InMemRecord { Solution sol; double objective; double time_sec; };
static vector<InMemRecord> g_history;

using clk = std::chrono::steady_clock;
inline double now_since(const clk::time_point &start) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(clk::now() - start).count();
}

void fillDist(const ProblemData& data) {
    int vcount = (int)data.villages.size();
    int ccount = (int)data.cities.size();
    dist_btw_villages.assign(vcount, vector<double>(vcount, 0.0));
    dist_btw_cities_villages.assign(ccount, vector<double>(vcount, 0.0));

    for (int i = 0; i < vcount; ++i)
        for (int j = 0; j < vcount; ++j)
            dist_btw_villages[i][j] = distance(data.villages[i].coords, data.villages[j].coords);

    for (int i = 0; i < ccount; ++i)
        for (int j = 0; j < vcount; ++j)
            dist_btw_cities_villages[i][j] = distance(data.cities[i], data.villages[j].coords);
}

double round_trip_distance(const vector<int>& village_indices, int city_idx) {
    if (village_indices.empty()) return 0.0;
    double d = dist_btw_cities_villages[city_idx][village_indices[0]];
    for (size_t k = 1; k < village_indices.size(); ++k)
        d += dist_btw_villages[village_indices[k-1]][village_indices[k]];
    d += dist_btw_cities_villages[city_idx][village_indices.back()];
    return d;
}

vector<int> find_reachable_villages_within_distance(int city_index, double radius) {
    vector<int> result;
    if (dist_btw_cities_villages.empty()) return result;
    int V = (int)dist_btw_cities_villages[0].size();
    for (int v = 0; v < V; ++v)
        if (dist_btw_cities_villages[city_index][v] <= radius)
            result.push_back(v);
    return result;
}

bool violates_caps(const Solution& sol, const ProblemData& problem) {
    int V = (int)problem.villages.size();
    vector<long long> meals(V, 0), other(V, 0);
    for (const auto& plan : sol) {
        for (const auto& trip : plan.trips) {
            for (const auto& d : trip.drops) {
                int vid = d.village_id - 1;
                if (vid < 0 || vid >= V) continue;
                meals[vid] += (long long)d.dry_food + (long long)d.perishable_food;
                other[vid] += (long long)d.other_supplies;
            }
        }
    }
    for (int i = 0; i < V; ++i) {
        long long meal_cap = 9LL * problem.villages[i].population;
        long long other_cap = 1LL * problem.villages[i].population;
        if (meals[i] > meal_cap) return true;
        if (other[i] > other_cap) return true;
    }
    return false;
}

bool violates_DMax(const Solution& sol, const ProblemData& problem) {
    for (const auto& plan : sol) {
        auto it = find_if(problem.helicopters.begin(), problem.helicopters.end(),
                          [&](const Helicopter& h){ return h.id == plan.helicopter_id; });
        if (it == problem.helicopters.end()) continue;
        const Helicopter& heli = *it;
        int city_index = heli.home_city_id - 1;
        double sumd = 0.0;
        for (const auto& trip : plan.trips) {
            vector<int> vids;
            for (const auto& d : trip.drops) {
                int idx = d.village_id - 1;
                if (idx >= 0 && idx < (int)problem.villages.size()) vids.push_back(idx);
            }
            sumd += round_trip_distance(vids, city_index);
        }
        if (sumd > problem.d_max + 1e-9) return true;
    }
    return false;
}

bool is_trip_feasible_tripwise(const Trip& trip, const Helicopter& heli, const ProblemData& problem) {
    double wsum = 0.0;
    for (const auto& d : trip.drops) {
        wsum += d.dry_food * problem.packages[0].weight
              + d.perishable_food * problem.packages[1].weight
              + d.other_supplies * problem.packages[2].weight;
    }
    if (wsum > heli.weight_capacity + 1e-9) return false;

    vector<int> vids;
    for (const auto& d : trip.drops) {
        int idx = d.village_id - 1;
        if (idx >= 0 && idx < (int)problem.villages.size()) vids.push_back(idx);
    }
    double rt = round_trip_distance(vids, heli.home_city_id - 1);
    if (rt > heli.distance_capacity + 1e-9) return false;
    return true;
}

bool violates_tripwise_constraints(const Solution& sol, const ProblemData& problem) {
    for (const auto& plan : sol) {
        auto it = find_if(problem.helicopters.begin(), problem.helicopters.end(),
                          [&](const Helicopter& h){ return h.id == plan.helicopter_id; });
        if (it == problem.helicopters.end()) continue;
        const Helicopter& heli = *it;
        for (const auto& trip : plan.trips) {
            if (!is_trip_feasible_tripwise(trip, heli, problem)) return true;
        }
    }
    return false;
}

double evaluate_solution(const Solution& solution, const ProblemData& problem) {
    struct Delivered { long long dry=0, wet=0, other=0; };
    vector<Delivered> delivered(problem.villages.size());
    unordered_map<int,double> heli_total_distance;
    double total_cost = 0.0;

    for (const auto& plan : solution) {
        auto it = find_if(problem.helicopters.begin(), problem.helicopters.end(),
                          [&](const Helicopter& h){ return h.id == plan.helicopter_id; });
        if (it == problem.helicopters.end()) continue;
        const Helicopter& heli = *it;

        double dist_sum = 0.0;
        for (const auto& trip : plan.trips) {
            vector<int> vidx;
            for (const auto& d : trip.drops) {
                int idx = d.village_id - 1;
                if (idx >= 0 && idx < (int)problem.villages.size()) {
                    delivered[idx].dry += d.dry_food;
                    delivered[idx].wet += d.perishable_food;
                    delivered[idx].other += d.other_supplies;
                    vidx.push_back(idx);
                }
            }
            dist_sum += round_trip_distance(vidx, heli.home_city_id - 1);
            total_cost += heli.fixed_cost;
        }
        total_cost += heli.alpha * dist_sum;
        heli_total_distance[heli.id] = dist_sum;
    }

    for (const auto& h : problem.helicopters) {
        double dsum = 0.0;
        auto it = heli_total_distance.find(h.id);
        if (it != heli_total_distance.end()) dsum = it->second;
        if (dsum > problem.d_max + 1e-9) return -1e15;
    }

    double total_value = 0.0;
    for (size_t i = 0; i < problem.villages.size(); ++i) {
        int pop = problem.villages[i].population;
        int max_meals = 9 * pop;
        long long wet_used = min<long long>(delivered[i].wet, max_meals);
        int rem = max_meals - (int)wet_used;
        long long dry_used = min<long long>(delivered[i].dry, rem);

        total_value += wet_used * problem.packages[1].value;
        total_value += dry_used * problem.packages[0].value;

        int other_cap = pop;
        long long other_used = min<long long>(delivered[i].other, other_cap);
        total_value += other_used * problem.packages[2].value;
    }

    return total_value - total_cost;
}

void trim_overshoot(Solution& sol, const ProblemData& problem) {
    int V = (int)problem.villages.size();
    vector<long long> delivered_meals(V, 0), delivered_other(V, 0);
    for (const auto& plan : sol) {
        for (const auto& trip : plan.trips) {
            for (const auto& d : trip.drops) {
                int vid = d.village_id - 1;
                if (vid < 0 || vid >= V) continue;
                delivered_meals[vid] += d.dry_food + d.perishable_food;
                delivered_other[vid] += d.other_supplies;
            }
        }
    }

    for (int v = 0; v < V; ++v) {
        long long meal_cap = 9LL * problem.villages[v].population;
        long long other_cap = 1LL * problem.villages[v].population;

        if (delivered_meals[v] > meal_cap) {
            long long excess = delivered_meals[v] - meal_cap;
            for (auto& plan : sol) {
                for (auto it_trip = plan.trips.rbegin(); it_trip != plan.trips.rend() && excess > 0; ++it_trip) {
                    for (auto it_drop = it_trip->drops.rbegin(); it_drop != it_trip->drops.rend() && excess > 0; ++it_drop) {
                        if (it_drop->village_id - 1 != v) continue;
                        int rem = min<long long>(excess, it_drop->dry_food);
                        it_drop->dry_food -= rem;
                        it_trip->dry_food_pickup -= rem;
                        excess -= rem;
                        delivered_meals[v] -= rem;
                        if (excess <= 0) break;
                        rem = min<long long>(excess, it_drop->perishable_food);
                        it_drop->perishable_food -= rem;
                        it_trip->perishable_food_pickup -= rem;
                        excess -= rem;
                        delivered_meals[v] -= rem;
                    }
                }
            }
        }

        if (delivered_other[v] > other_cap) {
            long long excess = delivered_other[v] - other_cap;
            for (auto& plan : sol) {
                for (auto it_trip = plan.trips.rbegin(); it_trip != plan.trips.rend() && excess > 0; ++it_trip) {
                    for (auto it_drop = it_trip->drops.rbegin(); it_drop != it_trip->drops.rend() && excess > 0; ++it_drop) {
                        if (it_drop->village_id - 1 != v) continue;
                        int rem = min<long long>(excess, it_drop->other_supplies);
                        it_drop->other_supplies -= rem;
                        it_trip->other_supplies_pickup -= rem;
                        excess -= rem;
                        delivered_other[v] -= rem;
                    }
                }
            }
        }
    }

    for (auto& plan : sol) {
        for (auto it = plan.trips.begin(); it != plan.trips.end();) {
            it->drops.erase(remove_if(it->drops.begin(), it->drops.end(),
                [](const Drop& d){ return d.dry_food==0 && d.perishable_food==0 && d.other_supplies==0; }),
                it->drops.end());
            if (it->drops.empty()) it = plan.trips.erase(it);
            else ++it;
        }
    }
}

static const double WET_PREF_P = 0.75;

Trip create_one_village_max_trip(const Helicopter& heli, const ProblemData& problem, mt19937& rng,
                                 vector<int>& remaining_meals, vector<int>& remaining_other,
                                 const clk::time_point &safe_deadline) {
    Trip best_trip{};
    int city_idx = heli.home_city_id - 1;

    struct PInfo { int idx; double weight; double value; double density; };
    vector<PInfo> pinfo(3);
    for (int i = 0; i < 3; ++i) {
        pinfo[i].idx = i;
        pinfo[i].weight = problem.packages[i].weight;
        pinfo[i].value  = problem.packages[i].value;
        pinfo[i].density = problem.packages[i].value / (problem.packages[i].weight + 1e-12);
    }
    sort(pinfo.begin(), pinfo.end(), [](const PInfo& a, const PInfo& b){ return a.density > b.density; });

    bernoulli_distribution use_wet_first(WET_PREF_P);
    double best_score = -1e18;
    int V = (int)problem.villages.size();

    for (int vidx = 0; vidx < V; ++vidx) {
        if (clk::now() >= safe_deadline) break;
        long long max_meals = remaining_meals[vidx];
        long long max_other = remaining_other[vidx];
        if (max_meals <= 0 && max_other <= 0) continue;

        const Village& v = problem.villages[vidx];
        double to_v = dist_btw_cities_villages[city_idx][vidx];
        double rt = 2.0 * to_v;
        if (rt > heli.distance_capacity + 1e-9) continue;

        long long allocate[3] = {0,0,0};
        double remaining_weight = heli.weight_capacity;
        bool wet_first = use_wet_first(rng);

        if (wet_first) {
            int pack_order[3] = {1,0,2};
            for (int p = 0; p < 3; ++p) {
                int type = pack_order[p];
                double w = problem.packages[type].weight;
                if (type == 2) {
                    long long can_by_weight = (long long)floor((remaining_weight + 1e-12) / w);
                    long long take = min(can_by_weight, max_other - allocate[2]);
                    if (take < 0) take = 0;
                    allocate[2] = take;
                    remaining_weight -= allocate[2] * w;
                } else {
                    long long allocated_meals_sofar = allocate[0] + allocate[1];
                    long long remaining_meals_allowed = max_meals - allocated_meals_sofar;
                    if (remaining_meals_allowed <= 0) continue;
                    long long can_by_weight = (long long)floor((remaining_weight + 1e-12) / w);
                    long long take = min(can_by_weight, remaining_meals_allowed);
                    if (take < 0) take = 0;
                    allocate[type] = take;
                    remaining_weight -= allocate[type] * w;
                }
            }
        } else {
            for (const auto &pi : pinfo) {
                int type = pi.idx;
                double w = pi.weight;
                if (type == 2) {
                    long long can_by_weight = (long long)floor((remaining_weight + 1e-12) / w);
                    long long take = min(can_by_weight, max_other - allocate[2]);
                    if (take < 0) take = 0;
                    allocate[2] = take;
                    remaining_weight -= allocate[2] * w;
                } else {
                    long long allocated_meals_sofar = allocate[0] + allocate[1];
                    long long remaining_meals_allowed = max_meals - allocated_meals_sofar;
                    if (remaining_meals_allowed <= 0) continue;
                    long long can_by_weight = (long long)floor((remaining_weight + 1e-12) / w);
                    long long take = min(can_by_weight, remaining_meals_allowed);
                    if (take < 0) take = 0;
                    allocate[type] = take;
                    remaining_weight -= allocate[type] * w;
                }
            }
        }

        if (allocate[0] + allocate[1] + allocate[2] == 0) continue;

        Trip cand{};
        cand.drops.push_back({v.id, (int)allocate[0], (int)allocate[1], (int)allocate[2]});
        cand.dry_food_pickup = (int)allocate[0];
        cand.perishable_food_pickup = (int)allocate[1];
        cand.other_supplies_pickup = (int)allocate[2];

        if (!is_trip_feasible_tripwise(cand, heli, problem)) {
            double w_total = allocate[0]*problem.packages[0].weight + allocate[1]*problem.packages[1].weight + allocate[2]*problem.packages[2].weight;
            if (w_total <= 0) continue;
            double scale = (heli.weight_capacity - 1e-9) / w_total;
            if (scale <= 0) continue;
            allocate[0] = (long long)floor(allocate[0] * scale);
            allocate[1] = (long long)floor(allocate[1] * scale);
            allocate[2] = (long long)floor(allocate[2] * scale);
            cand.drops.clear();
            cand.drops.push_back({v.id, (int)allocate[0], (int)allocate[1], (int)allocate[2]});
            cand.dry_food_pickup = (int)allocate[0];
            cand.perishable_food_pickup = (int)allocate[1];
            cand.other_supplies_pickup = (int)allocate[2];
            if (!is_trip_feasible_tripwise(cand, heli, problem)) continue;
        }

        double trip_value = allocate[0] * problem.packages[0].value
                          + allocate[1] * problem.packages[1].value
                          + allocate[2] * problem.packages[2].value;
        double trip_cost = heli.fixed_cost + heli.alpha * rt;
        double score = trip_value - trip_cost;

        if (score > best_score) {
            best_score = score;
            best_trip = cand;
        }
    }

    if (best_trip.drops.empty()) {
        Trip fallback{};
        int city_idx = heli.home_city_id - 1;
        auto reachable = find_reachable_villages_within_distance(city_idx, heli.distance_capacity);
        for (int vidx : reachable) {
            if (clk::now() >= safe_deadline) break;
            if (remaining_meals[vidx] <= 0 && remaining_other[vidx] <= 0) continue;
            const Village& v = problem.villages[vidx];
            int dry = min(remaining_meals[vidx], 1);
            int wet = 0;
            if (dry == 0 && remaining_other[vidx] > 0) {
                fallback.drops.push_back({v.id, 0, 0, 1});
                fallback.other_supplies_pickup = 1;
                break;
            }
            if (dry > 0) {
                fallback.drops.push_back({v.id, dry, wet, 0});
                fallback.dry_food_pickup = dry;
                break;
            }
        }
        if (is_trip_feasible_tripwise(fallback, heli, problem)) return fallback;
        return Trip{};
    }

    return best_trip;
}

void neighborhood_move(HelicopterPlan& plan, const ProblemData& problem, mt19937& rng,
                       const clk::time_point &safe_deadline) {
    if (clk::now() >= safe_deadline) return;
    if (plan.trips.empty()) return;

    auto plan_has_tripwise_violation = [&](const HelicopterPlan& p) -> bool {
        auto it = find_if(problem.helicopters.begin(), problem.helicopters.end(),
                          [&](const Helicopter& h){ return h.id == p.helicopter_id; });
        if (it == problem.helicopters.end()) return true;
        const Helicopter& heli = *it;
        for (const auto& t : p.trips) if (!is_trip_feasible_tripwise(t, heli, problem)) return true;
        return false;
    };

    uniform_int_distribution<> move_dist(0, 3);
    int mv = move_dist(rng);
    Trip& trip = plan.trips[0];

    switch (mv) {
        case 0: {
            int V = (int)problem.villages.size();
            if (V == 0) break;
            uniform_int_distribution<> pick_v(0, V-1);
            int chosen = pick_v(rng);
            uniform_int_distribution<> small(1,3);
            int r = small(rng);
            int add_dry=0, add_wet=0, add_other=0;
            if (r == 1) add_wet = 1;
            else if (r == 2) add_dry = 1;
            else add_other = 1;

            bool appended = false;
            for (auto &d : trip.drops) {
                if (d.village_id - 1 == chosen) {
                    d.dry_food += add_dry;
                    d.perishable_food += add_wet;
                    d.other_supplies += add_other;
                    trip.dry_food_pickup += add_dry;
                    trip.perishable_food_pickup += add_wet;
                    trip.other_supplies_pickup += add_other;
                    appended = true;
                    break;
                }
            }
            if (!appended) {
                trip.drops.push_back({problem.villages[chosen].id, add_dry, add_wet, add_other});
                trip.dry_food_pickup += add_dry;
                trip.perishable_food_pickup += add_wet;
                trip.other_supplies_pickup += add_other;
            }

            if (plan_has_tripwise_violation(plan)) {
                for (auto it = trip.drops.begin(); it != trip.drops.end(); ++it) {
                    if (it->village_id - 1 == chosen) {
                        it->dry_food -= add_dry;
                        it->perishable_food -= add_wet;
                        it->other_supplies -= add_other;
                        trip.dry_food_pickup -= add_dry;
                        trip.perishable_food_pickup -= add_wet;
                        trip.other_supplies_pickup -= add_other;
                        if (it->dry_food==0 && it->perishable_food==0 && it->other_supplies==0) trip.drops.erase(it);
                        break;
                    }
                }
            }
            break;
        }
        case 1: {
            if (trip.drops.size() <= 1) break;
            uniform_int_distribution<> pick(0, (int)trip.drops.size()-1);
            int idx = pick(rng);
            trip.dry_food_pickup -= trip.drops[idx].dry_food;
            trip.perishable_food_pickup -= trip.drops[idx].perishable_food;
            trip.other_supplies_pickup -= trip.drops[idx].other_supplies;
            trip.drops.erase(trip.drops.begin() + idx);
            break;
        }
        case 2: {
            if (trip.drops.empty()) break;
            uniform_int_distribution<> pick(0, (int)trip.drops.size()-1);
            int idx = pick(rng);
            Drop old = trip.drops[idx];
            Drop np = old;
            uniform_int_distribution<> adj(-3,3);
            np.dry_food = max(0, np.dry_food + adj(rng));
            np.perishable_food = max(0, np.perishable_food + adj(rng));
            np.other_supplies = max(0, np.other_supplies + adj(rng));
            trip.dry_food_pickup += (np.dry_food - old.dry_food);
            trip.perishable_food_pickup += (np.perishable_food - old.perishable_food);
            trip.other_supplies_pickup += (np.other_supplies - old.other_supplies);
            trip.drops[idx] = np;

            if (plan_has_tripwise_violation(plan)) {
                trip.drops[idx] = old;
                trip.dry_food_pickup -= (np.dry_food - old.dry_food);
                trip.perishable_food_pickup -= (np.perishable_food - old.perishable_food);
                trip.other_supplies_pickup -= (np.other_supplies - old.other_supplies);
            }
            break;
        }
        case 3: {
            if (trip.drops.size() < 2) break;
            shuffle(trip.drops.begin(), trip.drops.end(), rng);
            break;
        }
    }
}

void try_merge_pairs(HelicopterPlan& plan, const Helicopter& heli, const ProblemData& prob) {
    if (plan.trips.size() < 2) return;
    for (size_t i = 0; i < plan.trips.size(); ++i) {
        for (size_t j = i+1; j < plan.trips.size(); ++j) {
            Trip merged = plan.trips[i];
            for (const auto& d : plan.trips[j].drops) {
                merged.drops.push_back(d);
                merged.dry_food_pickup += d.dry_food;
                merged.perishable_food_pickup += d.perishable_food;
                merged.other_supplies_pickup += d.other_supplies;
            }
            if (is_trip_feasible_tripwise(merged, heli, prob)) {
                plan.trips[i] = merged;
                plan.trips.erase(plan.trips.begin() + j);
                return;
            }
        }
    }
}

void cleanup_merge(Solution& sol, const ProblemData& problem) {
    for (auto& plan : sol) {
        auto it = find_if(problem.helicopters.begin(), problem.helicopters.end(),
                          [&](const Helicopter& h){ return h.id == plan.helicopter_id; });
        if (it == problem.helicopters.end()) continue;
        bool changed = true;
        while (changed) {
            size_t before = plan.trips.size();
            try_merge_pairs(plan, *it, problem);
            changed = (plan.trips.size() < before);
        }
    }
}

Solution solve(const ProblemData& problem) {
    cerr << "solver started\n";
    fillDist(problem);

    random_device rd;
    mt19937 rng(rd());

    clk::time_point start = clk::now();
    double total_seconds = problem.time_limit_minutes * 60.0;
    long long total_ms = static_cast<long long>(total_seconds * 1000.0 + 0.5);
    clk::time_point hard_deadline = start + std::chrono::milliseconds(total_ms);
    clk::time_point safe_deadline = start + std::chrono::milliseconds((long long)(total_ms * 0.9));
    if (safe_deadline > hard_deadline) safe_deadline = hard_deadline;

    int V = (int)problem.villages.size();

    Solution best_solution;
    double best_score = -1e18;

    const int max_sideways = 8;

    long long per_restart_ms = std::max(100LL, total_ms / 10);

    for (int restart = 0; clk::now() < safe_deadline; ++restart) {
        if (clk::now() >= safe_deadline) break;
        clk::time_point restart_deadline = safe_deadline;

        vector<int> remaining_meals(V), remaining_other(V);
        for (int i = 0; i < V; ++i) {
            remaining_meals[i] = 9 * problem.villages[i].population;
            remaining_other[i] = problem.villages[i].population;
        }

        Solution current;
        for (const auto& h : problem.helicopters) {
            HelicopterPlan plan; plan.helicopter_id = h.id;
            Trip t = create_one_village_max_trip(h, problem, rng, remaining_meals, remaining_other, restart_deadline);
            if (!t.drops.empty()) {
                for (const auto& d : t.drops) {
                    int vid = d.village_id - 1;
                    if (vid >= 0 && vid < V) {
                        int meals = d.dry_food + d.perishable_food;
                        remaining_meals[vid] = max(0, remaining_meals[vid] - meals);
                        remaining_other[vid] = max(0, remaining_other[vid] - d.other_supplies);
                    }
                }
                plan.trips.push_back(t);
            }
            current.push_back(plan);
        }

        if (violates_tripwise_constraints(current, problem) || violates_caps(current, problem) || violates_DMax(current, problem)) {
            continue;
        }

        double current_score = evaluate_solution(current, problem);
        if (current_score > best_score + 1e-9) {
            best_score = current_score;
            best_solution = current;
            g_history.push_back({best_solution, best_score, now_since(start)});
            cerr << "best init restart " << restart << " score=" << best_score << " time=" << g_history.back().time_sec << "\n";
        }

        int sideways = 0;
        bool improved = true;
        while (improved) {
            if (clk::now() >= restart_deadline || clk::now() >= safe_deadline) { improved = false; break; }
            improved = false;

            for (size_t pidx = 0; pidx < current.size(); ++pidx) {
                if (clk::now() >= restart_deadline || clk::now() >= safe_deadline) { improved = false; break; }

                Solution trial = current;
                neighborhood_move(trial[pidx], problem, rng, restart_deadline);

                if (violates_tripwise_constraints(trial, problem) || violates_DMax(trial, problem) || violates_caps(trial, problem)) {
                    continue;
                }

                double trial_score = evaluate_solution(trial, problem);
                if (trial_score > current_score + 1e-9) {
                    current = std::move(trial);
                    current_score = trial_score;
                    improved = true;
                    sideways = 0;
                    if (current_score > best_score + 1e-9) {
                        best_score = current_score;
                        best_solution = current;
                        g_history.push_back({best_solution, best_score, now_since(start)});
                        cerr << "best found restart " << restart << " score=" << best_score << " time=" << g_history.back().time_sec << "\n";
                    }
                    break;
                } else if (fabs(trial_score - current_score) <= 1e-9 && sideways < max_sideways) {
                    current = std::move(trial);
                    current_score = trial_score;
                    improved = true;
                    sideways++;
                }
            }
        }

        cleanup_merge(current, problem);
        if (!violates_tripwise_constraints(current, problem) && !violates_DMax(current, problem) && !violates_caps(current, problem)) {
            current_score = evaluate_solution(current, problem);
            if (current_score > best_score + 1e-9) {
                best_score = current_score;
                best_solution = current;
                g_history.push_back({best_solution, best_score, now_since(start)});
                cerr << "best after cleanup " << restart << " score=" << best_score << " time=" << g_history.back().time_sec << "\n";
            }
        }

        if (clk::now() >= safe_deadline) break;
    }

    trim_overshoot(best_solution, problem);
    if (violates_tripwise_constraints(best_solution, problem) || violates_DMax(best_solution, problem) || violates_caps(best_solution, problem)) {
        trim_overshoot(best_solution, problem);
        cleanup_merge(best_solution, problem);
    }
    double final_score = evaluate_solution(best_solution, problem);
    g_history.push_back({best_solution, final_score, now_since(start)});
    cerr << "done. best = " << final_score << " history=" << g_history.size() << "\n";

    return best_solution;
}


