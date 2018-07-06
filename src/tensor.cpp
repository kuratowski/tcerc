#include "tensor.h"

#include <stdexcept>

namespace TCERC {

/* TODO: Change the implementation of tensor. */

Tensor::Tensor(const Subspace &subspace, int p, int q, const std::vector<iRRAM::REAL> &tensor) : subspace(subspace)
{
	// Any subspace validation here?

	// Any data validation here?
	int d = subspace.dimension();
	int expected_size = 1;
	for (int i = 0; i < p + q; i++) {
		expected_size *= d;
	}
	if (tensor.size() != expected_size) {
		throw std::invalid_argument("Input tensor size mismatch");
	}
	this->data = tensor;

	this->p = p;
	this->q = q;
	this->pq = std::vector<bool>();
	for (int i = 0; i < p; i++) {
		this->pq.push_back(false);
	}
	for (int i = 0; i < q; i++) {
		this->pq.push_back(true);
	}

	this->Amatrix = std::vector<std::vector<iRRAM::REAL>>();
	for (int i = 0; i < d; i++) {
		std::vector<iRRAM::REAL> row;
		for (int j = 0; j < d; j++) {
			row.push_back(i == j ? 1.0 : 0.0);
		}
		this->Amatrix.push_back(row);
	}

	this->subspace = subspace;
}

Tensor::Tensor(const Subspace &subspace, std::vector<bool> type, const std::vector<iRRAM::REAL> &tensor) : subspace(subspace)
{
	// Any subspace validation here?

	// Any data validation here?
	int d = subspace.dimension();
	int expected_size = 1;

	this->pq = type;
	this->p = 0;
	this->q = 0;
	for (int i = 0; i < type.size(); i++) {
		if (type[i])
			this->q++;
		else
			this->p++;
	}

	for (int i = 0; i < type.size(); i++) {
		expected_size *= d;
	}

	if (tensor.size() != expected_size) {
		throw std::invalid_argument("Input tensor size mismatch");
	}

	this->data = tensor;

	this->Amatrix = std::vector<std::vector<iRRAM::REAL>>();
	for (int i = 0; i < d; i++) {
		std::vector<iRRAM::REAL> row;
		for (int j = 0; j < d; j++) {
			row.push_back(i == j ? 1.0 : 0.0);
		}
		this->Amatrix.push_back(row);
	}
	this->subspace = subspace;

}

Tensor::Tensor(const Subspace &subspace, const iRRAM::REALMATRIX &matrix) : subspace(subspace)
{
	this->p = 1;
	this->q = 1;
	this->pq = {false, true};
	if (matrix.maxrow != subspace.dimension() || matrix.maxcolumn != subspace.dimension())
		throw std::invalid_argument("cannot create tensor with matrix");
	for (int i = 0; i < subspace.dimension(); i++) {
		for (int j = 0; j < subspace.dimension(); j++) {
			this->data.push_back(matrix(i, j));
		}
	}

	int d = subspace.dimension();
	this->Amatrix = std::vector<std::vector<iRRAM::REAL>>();
	for (int i = 0; i < d; i++) {
		std::vector<iRRAM::REAL> row;
		for (int j = 0; j < d; j++) {
			row.push_back(i == j ? 1.0 : 0.0);
		}
		this->Amatrix.push_back(row);
	}

	this->subspace = subspace;
}

bool Tensor::is_covariant(int k) const
{
	return this->pq[k];
}

const std::vector<bool> &Tensor::get_type() const
{
	return this->pq;
}

Tensor Tensor::operator+(const Tensor &t) const
{
	if (this->subspace != t.subspace || this->p != t.p || this->q != t.q || this->pq != t.pq)
		throw std::invalid_argument("Tensor of different subspace, or different order");

	std::vector<iRRAM::REAL> data;
	for (int i = 0; i < this->data.size(); i++) {
		data.push_back(this->data[i] + t.data[i]);
	}
	return Tensor(this->subspace, this->pq, data);
}

Tensor &Tensor::operator+=(const Tensor &t)
{
	if (this->subspace != t.subspace || this->p != t.p || this->q != t.q || this->pq != t.pq)
		throw std::invalid_argument("Tensor of different subspace, or different order");

	for (int i = 0; i < this->data.size(); i++) {
		this->data[i] += t.data[i];
	}

	return *this;
}

Tensor Tensor::operator-(const Tensor &t) const
{
	if (this->subspace != t.subspace || this->p != t.p || this->q != t.q || this->pq != t.pq)
		throw std::invalid_argument("Tensor of different subspace, or different order");

	std::vector<iRRAM::REAL> data;
	for (int i = 0; i < this->data.size(); i++) {
		data.push_back(this->data[i] - t.data[i]);
	}
	return Tensor(this->subspace, this->pq, data);
}

Tensor &Tensor::operator-=(const Tensor &t)
{
	if (this->subspace != t.subspace || this->p != t.p || this->q != t.q || this->pq != t.pq)
		throw std::invalid_argument("Tensor of different subspace, or different order");

	std::vector<iRRAM::REAL> data;
	for (int i = 0; i < this->data.size(); i++) {
		this->data[i] -= t.data[i];
	}

	return *this;
}

Tensor Tensor::operator*(const iRRAM::REAL &c) const
{
	std::vector<iRRAM::REAL> data;
	for (int i = 0; i < this->data.size(); i++) {
		data.push_back(this->data[i] * c);
	}
	return Tensor(this->subspace, this->pq, data);
}

Tensor &Tensor::operator*=(const iRRAM::REAL &c)
{
	for (int i = 0; i < this->data.size(); i++) {
		this->data[i] *= c;
	}
	return *this;
}

Tensor Tensor::operator/(const iRRAM::REAL &c) const
{
	std::vector<iRRAM::REAL> data;
	for (int i = 0; i < this->data.size(); i++) {
		data.push_back(this->data[i] / c);
	}
	return Tensor(this->subspace, this->pq, data);
}

Tensor &Tensor::operator/=(const iRRAM::REAL &c)
{
	for (int i = 0; i < this->data.size(); i++) {
		this->data[i] /= c;
	}
	return *this;
}

Tensor Tensor::operator*(const Tensor &t) const
{
	if (this->subspace != t.subspace)
		throw std::invalid_argument("Tensor of different subspace");
	int new_p = this->p + t.p;
	int new_q = this->q + t.q;
	std::vector<iRRAM::REAL> data;
	for (int i = 0; i < this->data.size(); i++) {
		for (int j = 0; j < t.data.size(); j++) {
			data.push_back(this->data[i] * t.data[j]);
		}
	}

	std::vector<bool> new_pq;
	for (int i = 0; i < this->pq.size(); i++) {
		new_pq.push_back(this->pq[i]);
	}
	for (int i = 0; i < t.pq.size(); i++) {
		new_pq.push_back(t.pq[i]);
	}

	return Tensor(this->subspace, new_pq, data);

}

Tensor &Tensor::contraction(int k, int l)
{

	int first = k > l ? l : k;
	int second = k > l ? k : l;
	if ((k < 0) || (k >= this->p + this->q) || 
			(l < 0) || (l >= this->p + this->q) || 
			(this->pq[first] == this->pq[second]))
		throw std::invalid_argument("Tensor cannot do contraction");

	int p = this->p;
	int q = this->q;
	int n = this->subspace.dimension();

	std::vector<iRRAM::REAL> data;

	std::vector<int> count;
	for (int i = 0; i < p + q - 2; i++)
		count.push_back(0);


	while (true) {

		iRRAM::REAL coeff = 0;
		for (int i = 0; i < n; i++) {

			std::vector<int> indices;
			for (int j = 0; j < p + q; j++){
				if (j == first) {
					indices.push_back(i);
				} else if (j == second) {
					indices.push_back(i);
				} else if (j < first) {
					indices.push_back(count[j]);
				} else if (j < second) {
					indices.push_back(count[j - 1]);
				} else {
					indices.push_back(count[j - 2]);
				}
			}
			coeff += this->coeff(indices);
		}
		data.push_back(coeff);

		int j;
		for (j = p + q - 3; j >= 0; j--) {
			count[j] ++;
			if (count[j] < n)
				break;
			count[j] = 0;
		}
		if (j < 0) {
			break;
		}

	}
	this->data = data;
	std::vector<bool> pq_contraction;
	for (int i = 0; i < p + q; i++) {
		if(i == first) continue;
		if(i == second) continue;
		pq_contraction.push_back(this->pq[i]);
	}
	this->pq = pq_contraction;
	this->p--;
	this->q--;
	return *this;
}

iRRAM::REAL Tensor::coeff(const std::vector<int> &indices) const 
{

	int p = this->p;
	int q = this->q;
	int n = this->subspace.dimension();

	if (indices.size() != (p + q)) 
		throw std::invalid_argument("Tensor invalid coefficient indices");

	int offset = 0;
	int step = 1;
	for (int i = 0; i < p + q; i++) {
		offset += step * indices[p + q - 1 - i];
		step *= n;
	}
		
	/*
	printf("Coeff : ");
	for (int i = 0; i < p + q; i++) {
		printf("[%d]", indices[i]);
	}
	printf(" -> ");
	iRRAM::cout << this->data[offset] << std::endl;
	*/

	return this->data[offset];

}

void Tensor::set_isomorphism(const std::vector<std::vector<iRRAM::REAL>> &Amatrix)
{
	this->Amatrix = Amatrix;
}

Tensor &Tensor::dualize(int index)
{
	// direction True => V->V*
	// direction False => V*->V
	if (this->p == 0)
		throw std::invalid_argument("Can't Dualize");
	std::vector<int> index_dualized;
	for (int i = 0; i < this->p + this->q; i++) {
		index_dualized.push_back(0);
	}

	std::vector<iRRAM::REAL> data;
	while (true) {
		iRRAM::REAL Amatrix_sum = 0;
		for (int c = 0; c < this->subspace.dimension(); c++) {
			std::vector<int> index_original_c;
			for (int i = 0; i < this->p + this->q; i++) {
				if (i == index) {
					index_original_c.push_back(c);
				} else {
					index_original_c.push_back(index_dualized[i]);
				}
			}
			int offset = 0;
			int power = 1;
			for (int i = 0; i < this->p + this->q; i++) {
				offset += index_original_c[this->p + this->q - 1 - i] * power;
				power *= this->subspace.dimension();
			}
			if (this->pq[index] == 0) { // V -> V*
				Amatrix_sum += this->data[offset] * this->Amatrix[c][index_dualized[index]];
			} else { // V* -> V
				Amatrix_sum += this->data[offset] * this->Amatrix[index_dualized[index]][c];
			}
		}
		data.push_back(Amatrix_sum);
		int i;
		for (i = this->p + this->q - 1; i >= 0; i--) {
			index_dualized[i] ++;
			if (index_dualized[i] == this->subspace.dimension()) {
				index_dualized[i] = 0;
			} else {
				break;
			}
		}
		if (i == -1) {
			break;
		}
	}
	this->data = data;
	if (this->pq[index] == 0) { // V -> V*
		this->p--;
		this->q++;
		this->pq[index] = true;
	} else { // V* -> V
		this->p++;
		this->q--;
		this->pq[index] = false;
	}
	return *this;
}

Tensor &Tensor::transpose(const std::vector<int> &permutation) 
{
	if (permutation.size() != this->p + this->q)
		throw std::invalid_argument("Invalid permutation");

	std::vector<int> permutation_inv;
	for (int i = 0; i < this->p + this->q; i++) {
		permutation_inv.push_back(-1);
	}
	for (int i = 0; i < this->p + this->q; i++) {
		permutation_inv[permutation[i]] = i;
	}

	std::vector<int> index_new;
	for (int i = 0; i < this->p + this->q; i++) {
		index_new.push_back(0);
	}

	std::vector<iRRAM::REAL> data;
	while (true) {
		// Apply inverse permutation to compute original index
		std::vector<int> index_original;
		for (int i = 0; i < this->p + this->q; i++) {
			index_original.push_back(-1);
		}
		for (int i = 0; i < this->p + this->q; i++) {
			index_original[permutation_inv[i]] = index_new[i];
		}

		// Convert original index vector to original offset
		int offset = 0;
		int powers = 1;
		for (int i = 0; i < this->p + this->q; i++) {
			offset += index_original[this->p + this->q - 1 - i] * powers;
		       	powers *= this->subspace.dimension();
		}

		data.push_back(this->data[offset]);
		
		int i;
		// Next for loop
		for (i = this->p + this->q - 1; i >= 0; i--) {
			index_new[i] ++;
			if (index_new[i] == this->subspace.dimension()) {
				index_new[i] = 0;
			} else {
				break;
			}
		}
		if (i == -1) {
			break;
		}

	}
	std::vector<bool> pq_permuted;
	for (int i = 0; i < this->p + this->q; i++) {
		pq_permuted.push_back(false);
	}
	for (int i = 0; i < this->p + this->q; i++) {
		pq_permuted[permutation[i]] = this->pq[i];
	}
	this->pq = pq_permuted;
	this->data = data;
	return *this;
}







} // TCERC
