class StaticMotionModel
{
	public:

		StaticMotionModel();

		inline gravis::d3Vector getPosChange(
			const std::vector<double>& x,
			int particle,
			int mode,
			int offset) const
		{
			return gravis::d3Vector(0.0, 0.0, 0.0);
		}
};
