################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../fichiers.cpp \
../fitness_diplo.cpp \
../main.cpp \
../migration_diplo.cpp \
../mutation_diplo.cpp \
../parameter.cpp \
../ranbin.cpp \
../rec.cpp \
../recursion_diplo.cpp 

OBJS += \
./fichiers.o \
./fitness_diplo.o \
./main.o \
./migration_diplo.o \
./mutation_diplo.o \
./parameter.o \
./ranbin.o \
./rec.o \
./recursion_diplo.o 

CPP_DEPS += \
./fichiers.d \
./fitness_diplo.d \
./main.d \
./migration_diplo.d \
./mutation_diplo.d \
./parameter.d \
./ranbin.d \
./rec.d \
./recursion_diplo.d 


# Each subdirectory must supply rules for building sources it contributes
fichiers.o: ../fichiers.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -Iconfig++ -O2 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"fichiers.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -Iconfig++ -O2 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


