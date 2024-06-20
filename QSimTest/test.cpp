#include "pch.h"

#include "..//QSim/Qbit.h"
#include "..//QSim/RM_RSA.h"
#include "..//QSim/Alg.h"

TEST(Qbit, normalization)
{
	Qbit<double> test(2);
	test[1] = 2;
	test[2] = 2;
	test.normalization();
	EXPECT_EQ(test[1], test[2]);
	EXPECT_FLOAT_EQ(norm(test[1]) + norm(test[2]), 1);
}

TEST(Qbit, attach)
{
	Qbit<double> test1(2);
	Qbit<double> test2(1);
	test1[0] = 1;
	test1[3] = 1;
	test2[0] = -1;
	test2[1] = 1;
	test1.normalization();
	test2.normalization();
	test1.attach(test2);
	Qbit<double> test3(4);
	test3[0] = -1;
	test3[3] = -1;
	test3[4] = 1;
	test3[7] = 1;
	test3.normalization();
	for (size_t i = 0; i < 8; i++)
	{
		EXPECT_TRUE(norm(test1[i] - test3[i]) < 0.00000001);
	}
}

TEST(Qbit, condition_exp)
{
	Qbit<double> test1(2);;
	test1[0] = 1;
	test1[3] = 1;
	test1.normalization();
	vector<size_t> res = test1.condition_exp(0, 2, 1000);
	EXPECT_TRUE(abs(int(res[0]) - int(res[3])) < 100);
	EXPECT_TRUE(res[0] + res[3] == 1000);
	res = test1.condition_exp(0, 2, 1000, true);
	EXPECT_TRUE(res[0] + res[3] == 1);
	if (res[0] == 1) {
		EXPECT_TRUE(test1[0] == 1.);
	} else{
		EXPECT_TRUE(test1[3] == 1.);
	}
}

TEST(Qbit_gate, gate_H)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test2[0] = 1;
	test2[1] = 1;
	test1.H(0);
	test2.normalization();
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit_gate, gate_X)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test2[1] = 1;
	test1.X(0);
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit_gate, gate_Y)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test2[1] = im;
	test1.Y(0);
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit_gate, gate_Z)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test1[1] = 1;
	test1.normalization();
	test2[0] = 1;
	test2[1] = -1;
	test2.normalization();
	test1.Z(0);
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit_gate, gate_P)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test1[1] = 1;
	test1.normalization();
	test2[0] = 1;
	test2[1] = im;
	test2.normalization();
	test1.P(PI / 4);
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_TRUE((norm(test1[i] - test2[i])) < 0.00000001);
	}
}

TEST(Qbit_gate, gate_CU)
{
	Qbit<double> test1(2);
	Qbit<double> test2(2);
	test1[2] = 1;
	test1[3] = 1;
	test1.normalization();
	test2[1] = 1;
	test2[2] = 1;
	test2.normalization();
	test1.CNOT(0, 1);
	for (size_t i = 0; i < 4; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit_gate, SWAP)
{
	Qbit<double> test1(3);
	Qbit<double> test2(3);
	test1[3] = 1;
	test1[4] = 1;
	test1.normalization();
	test2[1] = 1;
	test2[6] = 1;
	test2.normalization();
	test1.SWAP(0, 2);
	for (size_t i = 0; i < 8; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit_alg, QFT_RQFT)
{
	Qbit<double> test1(3);
	test1[0] = 1;
	test1.QFT(0, 3);
	test1.RQFT(0, 3);
	EXPECT_TRUE(1.0 - norm(test1[0]) < 0.00000001);
}

TEST(Qbit_alg, ADD_SUB)
{
	Qbit<double> test1(3);
	test1[3] = 1;
	test1.QFT(0, 3);
	test1.ADD(2, 0, 3);
	test1.RQFT(0, 3);
	EXPECT_TRUE(1.0 - norm(test1[5]) < 0.00000001);
	test1.QFT(0, 3);
	test1.SUB(3, 0, 3);
	test1.RQFT(0, 3);
	EXPECT_TRUE(1.0 - norm(test1[2]) < 0.00000001);
	test1.QFT(0, 3);
	test1.ADD(7, 0, 3);
	test1.RQFT(0, 3);
	EXPECT_TRUE(1.0 - norm(test1[1]) < 0.00000001);
	test1.QFT(1, 3);
	test1.SUB(2, 1, 3);
	test1.RQFT(1, 3);
	EXPECT_TRUE(1.0 - norm(test1[5]) < 0.00000001);
}

TEST(Qbit_alg, ADD_SUB_BY_MODULE)
{
	Qbit<double> test1(5);
	test1[1] = 1;
	test1.QFT(0, 4);
	test1.ModADD(1, 4, 0, 5);
	test1.RQFT(0, 4);
	EXPECT_TRUE(1.0 - norm(test1[2]) < 0.00000001);

	test1.QFT(0, 4);
	test1.ModADD(3, 4, 0, 5);
	test1.RQFT(0, 4);
	EXPECT_TRUE(1.0 - norm(test1[1]) < 0.00000001);

	test1.QFT(0, 4);
	test1.ModSUB(3, 4, 0, 5);
	test1.RQFT(0, 4);
	EXPECT_TRUE(1.0 - norm(test1[2]) < 0.00000001);
}

TEST(Qbit_alg, MUL_RMUL_MULX_BY_MODULE)
{
	Qbit<double> test1(5);
	test1[1] = 1;
	test1.QFT(0, 4);
	test1.ModMULX(2, 3, 4, 0, 5);
	test1.RQFT(0, 4);
	EXPECT_TRUE(1.0 - norm(test1[3]) < 0.00000001);

	Qbit<double> test2(6);
	test2[2] = 1;
	test2.QFT(2, 5);
	test2.ModMUL(2, 3, 0, 6);
	test2.RQFT(2, 5);
	EXPECT_TRUE(1.0 - norm(test2[2 + (1 << 2)]) < 0.00000001);

	Qbit<double> test3(6);
	test3[3] = 1;
	test3.QFT(2, 5);
	test3.RModMUL(3, 4, 0, 6);
	test3.RQFT(2, 5);
	EXPECT_TRUE(1.0 - norm(test3[3 + (3 << 2)]) < 0.00000001);
}

TEST(Qbit_alg, unitMUL)
{
	Qbit<double> test2(6);
	test2[3] = 1;
	test2.QFT(2, 5);
	test2.unitMUL(3, 4, 0, 6);
	test2.RQFT(2, 5);
	EXPECT_TRUE(1.0 - norm(test2[1]) < 0.00000001);
}

TEST(Qbit_alg, SHOR)
{
	Qbit<double> test2(18);
	test2.Shor(2, 15, 0, 18, 0);
	vector<size_t> temp = test2.condition_exp(0, 8, (1i64 << 11));
	for (size_t i = 0; i < 256; i++)
	{
		if (temp[i] != 0)
		{
			EXPECT_TRUE((i % 64) == 0);
		}
	}
}

TEST(RMR_SA, RSA)
{
	unsigned int m = 4;
	EXPECT_EQ(RSA_decryption(RSA_encryption(m)), m);
}

TEST(U_dec, H)
{
	Qbit<double> test1(1);
	test1[0] = 1;
	test1.SX(0);
	test1.RZ(0, PI/2);
	test1.SX(0);
	Qbit<double> test2(1);
	test2[0] = 1;
	test2[1] = 1;
	test2.normalization();
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_TRUE(norm(test1[i] - test2[i]) < 0.00000001);
	}
	test1.H(0);
	EXPECT_TRUE(norm(test1[0] - 1.) < 0.00000001);
}

TEST(U_dec, ZY)
{
	Qbit<double> test1(1);
	test1[0] = 1;
	test1.X(0);
	test1.RZ(0, PI);
	test1.SX(0);
	test1.RZ(0, -PI);
	test1.SX(0);
	test1.Y(0);
	test1.Z(0);
	EXPECT_TRUE(norm(test1[0] - 1.) < 0.00000001);
}

TEST(U_dec, RY)
{
	Qbit<double> test1(1);
	test1[0] = 1;
	for (size_t i = 1; i <= 12; i++)
	{
		test1.SX(0);
		test1.RZ(0, PI/i);
		test1.X(0);
		test1.SX(0);
		test1.RY(0, -PI / i);
		EXPECT_TRUE(norm(test1[0] - 1.) < 0.00000001);
	}
}

TEST(State_preparation, abs_ang)
{
	vector<complex<double>> vec{ im, 0, 0, -1 };
	vector<vector<double>> res{ {PI / 2}, {0, PI} };
	vector<vector<double>> res_ang = abs_ang(vec);
	
	for (size_t i = 0; i < res.size(); i++)
	{
		for (size_t j = 0; j < res[i].size(); j++)
		{
			EXPECT_DOUBLE_EQ(res[i][j], res_ang[i][j]);
		}
	}
}

TEST(State_preparation, phase_ang)
{
	vector<complex<double>> vec{ im, 0, 0, -1 };
	vector<double> res{ PI / 2, 0, 0, PI };
	vector<double> res_ang = phase_ang(vec);

	for (size_t i = 0; i < res.size(); i++)
	{
		EXPECT_DOUBLE_EQ(res[i], res_ang[i]);
	}
}

TEST(State_preparation, X_bit_mask)
{
	Qbit<double> cir(3);
	cir[0] = 1;
	X_bit_mask(cir, 5);
	EXPECT_TRUE(cir[5] == 1.);
}

TEST(State_preparation, State_preparation)
{
	vector<complex<double>> vec{ im, 0, 0, -1 };
	normalization_vec(vec);
	vector<complex<double>> res = State_preparation(vec);
	for (size_t i = 0; i < vec.size(); i++)
	{
		EXPECT_TRUE(norm(vec[i] - res[i]) < 0.00000001);
	}

	vector<complex<double>> vec2{ im, 0, 0, -1, 1, 2.*im, 3, 4 };
	normalization_vec(vec2);
	
	vector<complex<double>> res2 = State_preparation(vec2);
	for (size_t i = 0; i < vec2.size(); i++)
	{
		EXPECT_TRUE(norm(vec2[i] - res2[i]) < 0.00000001);
	}
}
