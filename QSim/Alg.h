#pragma once

size_t Shor(size_t a, size_t N, size_t n, double error = 0, bool print = false, string out_file = "Shor_result.txt")
{
    if (print)
    {
        cout << "P, CP, CPP error = " << error << "%\n";
    }

    auto start = std::chrono::steady_clock::now();

    Qbit<double> qb(4 * n + 2);

    qb.Shor(a, N, 0, 4 * n + 2, error);

    if (print)
    {
        qb.condition_exp_cout(0, 2 * n, (1i64 << (n * 2 + 3)));
    }

    size_t res = 1;

    //qb.cleaning_up_small_errors();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    if (print)
    {
        std::cout << "Shor " << a << ' ' << N << " time: " << elapsed_seconds.count() << "s\n";
    }

    ofstream out;
    out.open(out_file);

    if (out.is_open())
    {
        vector<size_t> temp = qb.condition_exp_in_file(0, 2 * n, (1i64 << (n * 2 + 3)), out);

        for (size_t i = 1; i < temp.size(); i++)
        {
            if ((temp[i] > (temp[0] / 2)) && (temp[i] > (temp[i + 1])) && (temp[i] > (temp[i - 1])))
            {
                res++;
            }
        }

        if (print)
        {
            cout << res << "\n\n";
        }

        out << a << ' ' << N << ' ' << error << ' ' << res << ' ' << round(elapsed_seconds.count()) << '\n';
    }
    else
    {
        vector<size_t> temp = qb.condition_exp(0, 2 * n, (1i64 << (n * 2 + 3)));

        for (size_t i = 1; i < temp.size(); i++)
        {
            if ((temp[i] > (temp[0] / 2)) && (temp[i] > (temp[i + 1])) && (temp[i] > (temp[i - 1])))
            {
                res++;
            }
        }

        if (print)
        {
            cout << res << "\n\n";
        }
    }

    out.close();

    return res;
}

size_t RSA_hacking(size_t e, size_t pq, size_t n1, size_t c)
{
    auto start = std::chrono::steady_clock::now();

    size_t res1 = 1;
    size_t e_ = e;

    while (true)
    {
        res1 = Shor(e, pq, n1, 0, false, "RSA_Shor_result.txt");

        size_t tempp = (gcd(MyPow(e, res1 / 2) + 1, pq));
        size_t tempm = (gcd(MyPow(e, res1 / 2) - 1, pq));
        if ((res1 & 1) || ((tempp == pq && tempp == 1) && (tempm == pq && tempm == 1)))
        {
            e++;
            while (gcd(e, pq) != 1)
            {
                e++;
                if (e >= pq)
                {
                    e = 2;
                }

            }

            if (e == e_)
            {
                cout << "Wrong pq\n";
                break;
            }

            res1 = 1;
        }
        else
        {
            break;
        }
    }

    size_t p = gcd(MyPow(e, res1 / 2) + 1, pq);
    if (p == pq || p == 1)
    {
        p = gcd(MyPow(e, res1 / 2) - 1, pq);
    }
    size_t q = pq / p;
    size_t fn = (p - 1) * (q - 1);

    size_t d = inverse_element_by_mod(e_, fn);

    size_t m_ = ModExp(c, d, pq);

    cout << m_ << '\n';

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "RSA hacking time: " << elapsed_seconds.count() << "s\n\n";

    return (m_);
}

void Dense_coding(Qbit<double> Dense, size_t count_exp = 100)//Dense 2,3 qbits = 0
{
    Dense.H(2);
    Dense.CNOT(2, 3);

    Dense.CNOT(1, 2);
    //Z
    Dense.H(2);
    Dense.CNOT(0, 2);
    Dense.H(2);

    Dense.CNOT(2, 3);
    Dense.H(2);

    Dense.condition_exp_cout(2, 4, count_exp);//Alice = Bob -> res = st
}

void Bells_inequality(Qbit<double> Bell, double a = 0, double b = 0, size_t count_exp = 100)//Bell(2 qubit) in 11
{
    Bell.H(1);
    Bell.CNOT(1, 0);

    Bell.RZ(0, a);
    Bell.RZ(1, b);

    Bell.H(0);
    Bell.H(1);

    for (int i = 0; i < 4; i++)
        cout << Bell[i] << '\n';

    cout << '\n';

    Bell.condition_exp_cout(0, 2, count_exp);
}

//for State preparation
template<typename T>
vector<vector<T>> abs_ang(vector<complex<T>> vec_)
{
    size_t step = vec_.size();
    vector<vector<T>> level_ang{};
    vector<T> vec{};

    for (size_t i = 0; i < vec_.size(); i++)
    {
        vec.push_back(abs(vec_[i]));
    }

    while (step != 1)
    {
        step = step >> 1;
        level_ang.push_back(vector<T>{});
        vector<T> temp_vec{};
        for (size_t j = 0; j < step; j++) {
            T a = vec[j];
            T b = vec[j + step];
            T c = sqrt(a * a + b * b);
            temp_vec.push_back(c);
            if (c != 0)
            {
                level_ang.back().push_back(2 * acos(a / c));
            }
            else {
                level_ang.back().push_back(0);
            }
            vec.clear();
            vec = temp_vec;
        }
    }

    reverse(level_ang.begin(), level_ang.end());

    return level_ang;
}

//for State preparation
template<typename T>
vector<T> phase_ang(vector<complex<T>> vec_)
{
    vector<T> level_ang{};
    
    for (size_t i = 0; i < vec_.size(); i++)
    {
        level_ang.push_back(arg(vec_[i]));
        if (level_ang.back() < 0) {
            level_ang.back() += 2 * PI;
        }
    }

    return level_ang;
}

template<typename T>
void X_bit_mask(Qbit<T>& cir, size_t x_mask, size_t start = 0) {
    while (x_mask != 0)
    {
        if (x_mask & 1)
        {
            cir.X(start);
        }
        start += 1;
        x_mask = x_mask >> 1;
    }
}

vector<complex<double>> State_preparation(vector<complex<double>> vec)
{
    size_t vec_size = vec.size();
    size_t N = log2(vec_size);
    Qbit<double> cir(N);
    cir[0] = 1;

    vector<vector<double>> level_angles = abs_ang(vec);
    vector<double> phase = phase_ang(vec);

    vector<double> angles = level_angles[0];

    if (angles[0] != 0) {
        cir.RY(0, angles[0]);
    }
    if (N == 1) {
        if (phase[1] != 0)
        {
            cir.P(0 ,phase[1]);
        }
        if (phase[0] != 0)
        {
            cir.X(0);
            cir.P(0, phase[0]);
            cir.X(0);
        }
    }
    vector<size_t> temp{ 0 };
    
    for (size_t i = 1; i < N; i++)
    {
        size_t p_2 = (1 << i);
        angles.clear();
        angles = level_angles[i];
        size_t last_0 = p_2 - 1;
        for (size_t j = 0; j < p_2; j++)
        {  
            
            if (angles[j] != 0 || i == N - 1 || j == 0)
            {
                X_bit_mask(cir, ((p_2 - j - 1) ^ (p_2 - last_0 - 1)));
                last_0 = j;
            }

            if (i == N - 1 && phase[j] != 0)
            {
                cir.RZ(i, -phase[j]);
            }
            if (i == N - 1 && phase[j + p_2] != 0)
            {
                cir.X(i);
                cir.RZ(i, phase[j + p_2]);
                cir.X(i);
            }

            cir.RY(i, angles[j] / 2);
            cir.CnNOT(temp, i);
            cir.RY(i, -angles[j] / 2);

            if (i == N - 1 && phase[j] != 0)
            {
                cir.RZ(i, phase[j]);
            }

            cir.CnNOT(temp, i);

            if (i == N - 1 && phase[j + p_2] != 0)
            {
                cir.X(i);
                cir.RZ(i, -phase[j + p_2]);
                cir.X(i);
            }
        }

        X_bit_mask(cir, p_2 - last_0 - 1);

        temp.push_back(i);
    }

    return cir.get_data();
}
