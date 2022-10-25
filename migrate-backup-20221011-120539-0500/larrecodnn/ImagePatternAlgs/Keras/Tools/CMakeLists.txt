include(lar::PointIdAlgorithm)

cet_build_plugin(PointIdAlgKeras lar::PointIdAlgorithm
  LIBRARIES PRIVATE
  larrecodnn::ImagePatternAlgs_Keras
)

install_source()
