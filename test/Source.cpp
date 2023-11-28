#include <iostream>
#include <vector>

using namespace std;

int main()
{
	std::vector<int> a;
	a.resize(5);
	a = { 1, 2, 3, 4, 5 };
	cout << a[2] << endl;



	return 0;
}