#ifdef TORCH_ENABLED
#include <torch/torch.h>
#endif //TORCH_ENABLED

#include <iostream>

int main(int argc, char **argv)
{
#ifdef TORCH_ENABLED
		torch::Tensor tensor = torch::rand({2, 3});
		std::cout << tensor << std::endl;
#endif //TORCH_ENABLED

	return 0;
}

