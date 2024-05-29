PROJECT(FEM)
SET ( LIB_NAME		${PROJECT_NAME}   	)

#== Necessário para incluir arquivos .mod ======
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/fontes	)

#SET(PKG_COMMON_SRC 
SET(PKG_FEM_SRC 
	estruturaCSR.F90
	algMatricial.F90
	pardisoSistema.F90
	HYPRESistema.F90
	malha.F90
	utilitarios.F90
	variaveisGlobais.F90
	mInputReader.F90
	fluxo.F90
	leituraEscrita.F90
	funcoesDeForma.F90
	potencial.F90	
)


SET(EXT_LIBS
#	vtkCommon	
)

ADD_LIBRARY(${PROJECT_NAME} ${PKG_FEM_SRC})

FOREACH(c ${EXT_LIBS})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${c})
ENDFOREACH(c)

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}
#    DEBUG_POSTFIX "_${VS_DIR}_${PLATFORM}_d"
#    RELEASE_POSTFIX "_${VS_DIR}_${PLATFORM}
 )

