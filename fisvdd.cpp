#include <iostream>
#include <string>
#include <vector>
#include <glob.h>

#include <typeinfo> // for check variable type
#include <cmath> // sqrt
#include <numeric> // inner product
#include <tuple> // multiple return
#include <algorithm> // min max

using namespace std;


void print_2d_vector(vector<vector<double>> A){
  cout << "[";
  for (int i = 0; i < A.size(); i++){
    cout << "[";
    for (int j = 0; j < A[0].size(); j++){
      printf("%.5f ", A[i][j]);
    }
    cout << "]" << endl;
  }
}

void print_1d_vector(vector<double> A){
  cout << "[";
  for (int i = 0; i < A.size(); i++){
    printf("%.5f ", A[i]);
  }
  cout << "]" << endl;
}

double dot_product(vector<double> a, vector<double> b){
  double dot = 0;
  for (int i = 0; i < a.size(); i++)
    dot += a[i] * b[i];
  return dot;
}

vector<double> dot_product(vector<vector<double>> a, vector<double> b){
  vector<double> dot;
  for (int i = 0; i < a.size(); i++)
    dot.push_back(dot_product(a[i], b));
  return dot;
}

vector<vector<double>> outer_product(vector<double> a, vector<double> b){
  vector<vector<double>> outer;
  for (int i = 0; i < a.size(); i++){
    vector<double> outer_row;
    for (int j = 0; j < b.size(); j++){
      outer_row.push_back(a[i] * b[j]);
    }
    outer.push_back(outer_row);
  }
  return outer;
}

vector<double> concat_1d_vector(vector<vector<double>> vector_list){
  vector<double> concated;
  for (auto& n: vector_list) 
    concated.insert(end(concated), begin(n), end(n));
  return concated;
}

vector<vector<double>> concat_2d_vector(vector<vector<vector<double>>> vector_list){
  vector<vector<double>> concated;
  for (auto& n: vector_list) 
    concated.insert(end(concated), begin(n), end(n));
  return concated;
}

class fisvdd{
  public:
    double eps_cp;
    double eps_ol;
    vector<vector<double>> inv_A;
    vector<double> alpha;
    vector<vector<double>> sv;
    vector<double> obj_val;
    vector<double> init_sv;
    double score;
    double sigma;
    bool verbose = false;

    fisvdd(vector<double> _init_sv, double _sigma){
      sigma = _sigma;
      eps_cp = 1e-8;
      eps_ol = 1e-8;

      vector<double> init(1, 1);

      inv_A.push_back(init);
      alpha = init;
      sv.push_back(_init_sv);
      score = 1;
    }

    //void find_sv(vector<vector<double>> data);
    tuple<vector<double>, vector<vector<double>>> find_sv(vector<vector<double>> data);
    tuple<double, vector<double>> score_fcn(vector<double> new_data);
    vector<double> kernel(vector<double> new_data, vector<vector<double>> support_vector_sets, double sigma);
    void expand(vector<double> new_sv, vector<double> new_sim_vec);
    vector<vector<double>> up_inv(vector<vector<double>> prev_inv, vector<double> v);
    vector<vector<double>> shrink(void);
    vector<vector<double>> perm(vector<vector<double>> A,  int ind);
    vector<vector<double>> down_inv(vector<vector<double>> next_inv);
    void model_update(void);
};

tuple<vector<double>, vector<vector<double>>> fisvdd::find_sv(vector<vector<double>> data){
  vector<double> new_data;
  tuple<double, vector<double>> score_ret;
  vector<double> sim_vec; double new_score;

  for (int i = 1; i < data.size(); i++){
    new_data = data[i];
    score_ret = score_fcn(new_data);
    new_score = get<0>(score_ret); sim_vec = get<1>(score_ret);
    if (new_score > 0){
      expand(new_data, sim_vec);

      if (*min_element(alpha.begin(), alpha.end()) < 0){
        vector<vector<double>> backup = shrink();
        for (auto each: backup){
          score_ret = score_fcn(each);
          new_score = get<0>(score_ret); sim_vec = get<1>(score_ret);
          if (new_score > 0)
            expand(each, sim_vec);
        }
      }
      model_update();
    }
    if (verbose){
      printf("new_data:\n");
      print_1d_vector(new_data);
      printf("score: %.2f\n", new_score);
      printf("alpha(%d): ", i-1);
      print_1d_vector(alpha);
      cout << endl;
      if (i > 100 + 1) break;
    }
  }
  return tuple<vector<double>, vector<vector<double>>> (alpha, sv);
}


vector<vector<double>> fisvdd::up_inv(vector<vector<double>> prev_inv, vector<double> v){
  vector<double> p;
  vector<vector<double>> A;
  vector<double> C;
  double beta;
  p = dot_product(prev_inv, v);
  beta = 1 - dot_product(v, p);

  A = outer_product(p, p);
  for (int i = 0; i < A.size(); i++)
    for (int j = 0; j < A[0].size(); j++)
      A[i][j] = prev_inv[i][j] + A[i][j] / beta;

  for (int i = 0; i < p.size(); i++){
    A[i].push_back(-p[i] / beta); // np.hstack((A, B))
    C.push_back(-p[i] / beta);
  }
  C.push_back(1 / beta); // np.hstack((C, D))
  A.push_back(C); // np.vstack(())

  if (verbose){
    cout << "beta: "<< beta << endl;
    printf("up_inv_A Size: (%lo, %lo)\n",A.size(), A[0].size());
    print_2d_vector(A);
    printf("C[0]: %f, C Size: (%lo)\n", C[0], C.size());
  }

  return A; //res
}

vector<vector<double>> fisvdd::down_inv(vector<vector<double>> next_inv){
  int n = next_inv.size();
  double lamb = next_inv[n-1][n-1];
  if (verbose)
    printf("lamb: %f\n", lamb);
  vector<double> u = next_inv[n-1];
  u.pop_back();
  vector<vector<double>> u_outer = outer_product(u, u);
  vector<vector<double>> res(n - 1, vector<double> (n - 1));
  for (int i = 0; i < next_inv.size() - 1; i++)
    for (int j = 0; j < next_inv.size() - 1; j++)
      res[i][j] = next_inv[i][j] - u_outer[i][j] / lamb;

  return res;
}

void fisvdd::expand(vector<double> new_sv, vector<double> new_sim_vec){
  inv_A = up_inv(inv_A, new_sim_vec);
  alpha.clear();
  for (int i = 0; i < inv_A.size(); i++){
    double tmp_alpha = 0;
    for (int j = 0; j < inv_A[0].size(); j++){
      tmp_alpha += inv_A[i][j];
    }
    alpha.push_back(tmp_alpha);
  }
  sv.push_back(new_sv);

  if (verbose){
    cout << "alpha: " ;
    print_1d_vector(alpha);
  }
}

vector<vector<double>> fisvdd::shrink(void){
  vector<vector<double>> backup;
  while (true){
    int min_ind = min_element(alpha.begin(), alpha.end()) - alpha.begin();
    double min_val = *min_element(alpha.begin(), alpha.end());
    vector<double> data_out = sv[min_ind];
    backup.push_back(data_out);
    for (int i = 0; i < alpha.size(); i++)
      if (alpha[i] <= min_val){
        sv.erase(sv.begin() + i);
      }
    inv_A = perm(inv_A, min_ind);
    inv_A = down_inv(inv_A);
    alpha.clear();
    for (int i = 0; i < inv_A.size(); i++)
      alpha.push_back(accumulate(inv_A[i].begin(), inv_A[i].end(), 0.0));

    if (verbose){
      printf("minInd in shrink: %d\n", min_ind);
      printf("Alpha after shrink\n");
      print_1d_vector(alpha);
      printf("sv After shrink\n");
      print_2d_vector(sv);
    }
    if (*min_element(alpha.begin(), alpha.end()) > 0)
      break;
  }
  return backup;
}

vector<vector<double>> fisvdd::perm(vector<vector<double>> A, int ind){
  /* Permutation function */
  int n = A[0].size();
  vector<double> tmp = A[ind];
  A.erase(A.begin() + ind);
  A.push_back(tmp);
  for (int i = 0; i < A.size(); i++){
    double tmp = A[i][ind];
    A[i].erase(A[i].begin() + ind);
    A[i].push_back(tmp);
  }
  return A;
}


tuple<double, vector<double>> fisvdd::score_fcn(vector<double> new_data){
  vector<double> cur_sim_vec;
  double m;
  double res;
  cur_sim_vec = kernel(new_data, sv, sigma);
  if (verbose){
    printf("previous sv:\n");
    print_2d_vector(sv);

    printf("cur_sim_vec:\n");
    print_1d_vector(cur_sim_vec);
  }

  m = *max_element(cur_sim_vec.begin(), cur_sim_vec.end());
  if (m < eps_ol || m > 1 - eps_cp)
    res = -1;
  else
    res = score - dot_product(alpha, cur_sim_vec);

  return tuple<double, vector<double>>(res, cur_sim_vec);
}


vector<double> fisvdd::kernel(vector<double> new_data, vector<vector<double>> support_vector_sets, double sigma){
  vector<double> K;
  for(int i=0; i < int(support_vector_sets.size()); i++){
    double dist_sq = 0;
    for(int j = 0; j < int(support_vector_sets[i].size()); j++){
      double tmp = new_data[j] - support_vector_sets[i][j];
      dist_sq += tmp * tmp;
    }
    K.push_back(exp(-dist_sq / (2.0 * sigma * sigma)));
  }
  return K;
}

void fisvdd::model_update(void){
  double alpha_sum = 0.;
  for (auto& n: alpha) alpha_sum += n;
  score = 1 / alpha_sum;
  for (int i = 0; i < alpha.size(); i++)
    alpha[i] = alpha[i] / alpha_sum;
}


void print_feature(vector<vector<double>> train_feats){
  cout << train_feats.size() << ' ' << train_feats[0].size() << endl;
  for (int i = 0; i < int(train_feats[0].size()); i++){
    cout << train_feats[0][i] << endl; 
  }
}

