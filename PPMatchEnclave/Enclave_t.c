#include "Enclave_t.h"

#include "sgx_trts.h" /* for sgx_ocalloc, sgx_is_outside_enclave */
#include "sgx_lfence.h" /* for sgx_lfence */

#include <errno.h>
#include <mbusafecrt.h> /* for memcpy_s etc */
#include <stdlib.h> /* for malloc/free etc */

#define CHECK_REF_POINTER(ptr, siz) do {	\
	if (!(ptr) || ! sgx_is_outside_enclave((ptr), (siz)))	\
		return SGX_ERROR_INVALID_PARAMETER;\
} while (0)

#define CHECK_UNIQUE_POINTER(ptr, siz) do {	\
	if ((ptr) && ! sgx_is_outside_enclave((ptr), (siz)))	\
		return SGX_ERROR_INVALID_PARAMETER;\
} while (0)

#define CHECK_ENCLAVE_POINTER(ptr, siz) do {	\
	if ((ptr) && ! sgx_is_within_enclave((ptr), (siz)))	\
		return SGX_ERROR_INVALID_PARAMETER;\
} while (0)

#define ADD_ASSIGN_OVERFLOW(a, b) (	\
	((a) += (b)) < (b)	\
)


typedef struct ms_ecall_init_Query_t {
	int* ms_Query;
	size_t ms_matrixlen;
	int* ms_Label;
	size_t ms_labellen;
	int ms_querysize;
	int ms_labelsize;
} ms_ecall_init_Query_t;

typedef struct ms_ecall_BF_t {
	char* ms_bits;
	size_t ms_bitlen;
	char* ms_result;
	size_t ms_resultlen;
	int ms_numHashs;
	int ms_querysize;
	int ms_centerlabel;
} ms_ecall_BF_t;

typedef struct ms_ocall_print_string_t {
	const char* ms_str;
} ms_ocall_print_string_t;

typedef struct ms_sgx_oc_cpuidex_t {
	int* ms_cpuinfo;
	int ms_leaf;
	int ms_subleaf;
} ms_sgx_oc_cpuidex_t;

typedef struct ms_sgx_thread_wait_untrusted_event_ocall_t {
	int ms_retval;
	const void* ms_self;
} ms_sgx_thread_wait_untrusted_event_ocall_t;

typedef struct ms_sgx_thread_set_untrusted_event_ocall_t {
	int ms_retval;
	const void* ms_waiter;
} ms_sgx_thread_set_untrusted_event_ocall_t;

typedef struct ms_sgx_thread_setwait_untrusted_events_ocall_t {
	int ms_retval;
	const void* ms_waiter;
	const void* ms_self;
} ms_sgx_thread_setwait_untrusted_events_ocall_t;

typedef struct ms_sgx_thread_set_multiple_untrusted_events_ocall_t {
	int ms_retval;
	const void** ms_waiters;
	size_t ms_total;
} ms_sgx_thread_set_multiple_untrusted_events_ocall_t;

static sgx_status_t SGX_CDECL sgx_ecall_lambdas_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_lambdas_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_init_Query(void* pms)
{
	CHECK_REF_POINTER(pms, sizeof(ms_ecall_init_Query_t));
	//
	// fence after pointer checks
	//
	sgx_lfence();
	ms_ecall_init_Query_t* ms = SGX_CAST(ms_ecall_init_Query_t*, pms);
	sgx_status_t status = SGX_SUCCESS;
	int* _tmp_Query = ms->ms_Query;
	size_t _tmp_matrixlen = ms->ms_matrixlen;
	size_t _len_Query = _tmp_matrixlen;
	int* _in_Query = NULL;
	int* _tmp_Label = ms->ms_Label;
	size_t _tmp_labellen = ms->ms_labellen;
	size_t _len_Label = _tmp_labellen;
	int* _in_Label = NULL;

	CHECK_UNIQUE_POINTER(_tmp_Query, _len_Query);
	CHECK_UNIQUE_POINTER(_tmp_Label, _len_Label);

	//
	// fence after pointer checks
	//
	sgx_lfence();

	if (_tmp_Query != NULL && _len_Query != 0) {
		if ( _len_Query % sizeof(*_tmp_Query) != 0)
		{
			status = SGX_ERROR_INVALID_PARAMETER;
			goto err;
		}
		_in_Query = (int*)malloc(_len_Query);
		if (_in_Query == NULL) {
			status = SGX_ERROR_OUT_OF_MEMORY;
			goto err;
		}

		if (memcpy_s(_in_Query, _len_Query, _tmp_Query, _len_Query)) {
			status = SGX_ERROR_UNEXPECTED;
			goto err;
		}

	}
	if (_tmp_Label != NULL && _len_Label != 0) {
		if ( _len_Label % sizeof(*_tmp_Label) != 0)
		{
			status = SGX_ERROR_INVALID_PARAMETER;
			goto err;
		}
		_in_Label = (int*)malloc(_len_Label);
		if (_in_Label == NULL) {
			status = SGX_ERROR_OUT_OF_MEMORY;
			goto err;
		}

		if (memcpy_s(_in_Label, _len_Label, _tmp_Label, _len_Label)) {
			status = SGX_ERROR_UNEXPECTED;
			goto err;
		}

	}

	ecall_init_Query(_in_Query, _tmp_matrixlen, _in_Label, _tmp_labellen, ms->ms_querysize, ms->ms_labelsize);

err:
	if (_in_Query) free(_in_Query);
	if (_in_Label) free(_in_Label);
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_BF(void* pms)
{
	CHECK_REF_POINTER(pms, sizeof(ms_ecall_BF_t));
	//
	// fence after pointer checks
	//
	sgx_lfence();
	ms_ecall_BF_t* ms = SGX_CAST(ms_ecall_BF_t*, pms);
	sgx_status_t status = SGX_SUCCESS;
	char* _tmp_bits = ms->ms_bits;
	size_t _tmp_bitlen = ms->ms_bitlen;
	size_t _len_bits = _tmp_bitlen;
	char* _in_bits = NULL;
	char* _tmp_result = ms->ms_result;
	size_t _tmp_resultlen = ms->ms_resultlen;
	size_t _len_result = _tmp_resultlen;
	char* _in_result = NULL;

	CHECK_UNIQUE_POINTER(_tmp_bits, _len_bits);
	CHECK_UNIQUE_POINTER(_tmp_result, _len_result);

	//
	// fence after pointer checks
	//
	sgx_lfence();

	if (_tmp_bits != NULL && _len_bits != 0) {
		if ( _len_bits % sizeof(*_tmp_bits) != 0)
		{
			status = SGX_ERROR_INVALID_PARAMETER;
			goto err;
		}
		_in_bits = (char*)malloc(_len_bits);
		if (_in_bits == NULL) {
			status = SGX_ERROR_OUT_OF_MEMORY;
			goto err;
		}

		if (memcpy_s(_in_bits, _len_bits, _tmp_bits, _len_bits)) {
			status = SGX_ERROR_UNEXPECTED;
			goto err;
		}

	}
	if (_tmp_result != NULL && _len_result != 0) {
		if ( _len_result % sizeof(*_tmp_result) != 0)
		{
			status = SGX_ERROR_INVALID_PARAMETER;
			goto err;
		}
		if ((_in_result = (char*)malloc(_len_result)) == NULL) {
			status = SGX_ERROR_OUT_OF_MEMORY;
			goto err;
		}

		memset((void*)_in_result, 0, _len_result);
	}

	ecall_BF(_in_bits, _tmp_bitlen, _in_result, _tmp_resultlen, ms->ms_numHashs, ms->ms_querysize, ms->ms_centerlabel);
	if (_in_result) {
		if (memcpy_s(_tmp_result, _len_result, _in_result, _len_result)) {
			status = SGX_ERROR_UNEXPECTED;
			goto err;
		}
	}

err:
	if (_in_bits) free(_in_bits);
	if (_in_result) free(_in_result);
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_auto_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_auto_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_decltype_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_decltype_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_strongly_typed_enum_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_strongly_typed_enum_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_range_based_for_loops_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_range_based_for_loops_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_static_assert_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_static_assert_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_virtual_function_control_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_virtual_function_control_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_delegating_constructors_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_delegating_constructors_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_std_function_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_std_function_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_cxx11_algorithms_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_cxx11_algorithms_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_variadic_templates_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_variadic_templates_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_SFINAE_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_SFINAE_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_initializer_list_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_initializer_list_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_rvalue_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_rvalue_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_nullptr_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_nullptr_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_enum_class_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_enum_class_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_new_container_classes_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_new_container_classes_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_tuple_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_tuple_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_shared_ptr_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_shared_ptr_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_atomic_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_atomic_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_mutex_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_mutex_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_print_final_value_mutex_demo(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_print_final_value_mutex_demo();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_mutex_demo_no_protection(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_mutex_demo_no_protection();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_print_final_value_no_protection(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_print_final_value_no_protection();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_condition_variable_run(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_condition_variable_run();
	return status;
}

static sgx_status_t SGX_CDECL sgx_ecall_condition_variable_load(void* pms)
{
	sgx_status_t status = SGX_SUCCESS;
	if (pms != NULL) return SGX_ERROR_INVALID_PARAMETER;
	ecall_condition_variable_load();
	return status;
}

SGX_EXTERNC const struct {
	size_t nr_ecall;
	struct {void* ecall_addr; uint8_t is_priv; uint8_t is_switchless;} ecall_table[28];
} g_ecall_table = {
	28,
	{
		{(void*)(uintptr_t)sgx_ecall_lambdas_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_init_Query, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_BF, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_auto_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_decltype_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_strongly_typed_enum_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_range_based_for_loops_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_static_assert_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_virtual_function_control_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_delegating_constructors_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_std_function_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_cxx11_algorithms_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_variadic_templates_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_SFINAE_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_initializer_list_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_rvalue_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_nullptr_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_enum_class_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_new_container_classes_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_tuple_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_shared_ptr_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_atomic_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_mutex_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_print_final_value_mutex_demo, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_mutex_demo_no_protection, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_print_final_value_no_protection, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_condition_variable_run, 0, 0},
		{(void*)(uintptr_t)sgx_ecall_condition_variable_load, 0, 0},
	}
};

SGX_EXTERNC const struct {
	size_t nr_ocall;
	uint8_t entry_table[6][28];
} g_dyn_entry_table = {
	6,
	{
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
	}
};


sgx_status_t SGX_CDECL ocall_print_string(const char* str)
{
	sgx_status_t status = SGX_SUCCESS;
	size_t _len_str = str ? strlen(str) + 1 : 0;

	ms_ocall_print_string_t* ms = NULL;
	size_t ocalloc_size = sizeof(ms_ocall_print_string_t);
	void *__tmp = NULL;


	CHECK_ENCLAVE_POINTER(str, _len_str);

	if (ADD_ASSIGN_OVERFLOW(ocalloc_size, (str != NULL) ? _len_str : 0))
		return SGX_ERROR_INVALID_PARAMETER;

	__tmp = sgx_ocalloc(ocalloc_size);
	if (__tmp == NULL) {
		sgx_ocfree();
		return SGX_ERROR_UNEXPECTED;
	}
	ms = (ms_ocall_print_string_t*)__tmp;
	__tmp = (void *)((size_t)__tmp + sizeof(ms_ocall_print_string_t));
	ocalloc_size -= sizeof(ms_ocall_print_string_t);

	if (str != NULL) {
		ms->ms_str = (const char*)__tmp;
		if (_len_str % sizeof(*str) != 0) {
			sgx_ocfree();
			return SGX_ERROR_INVALID_PARAMETER;
		}
		if (memcpy_s(__tmp, ocalloc_size, str, _len_str)) {
			sgx_ocfree();
			return SGX_ERROR_UNEXPECTED;
		}
		__tmp = (void *)((size_t)__tmp + _len_str);
		ocalloc_size -= _len_str;
	} else {
		ms->ms_str = NULL;
	}
	
	status = sgx_ocall(0, ms);

	if (status == SGX_SUCCESS) {
	}
	sgx_ocfree();
	return status;
}

sgx_status_t SGX_CDECL sgx_oc_cpuidex(int cpuinfo[4], int leaf, int subleaf)
{
	sgx_status_t status = SGX_SUCCESS;
	size_t _len_cpuinfo = 4 * sizeof(int);

	ms_sgx_oc_cpuidex_t* ms = NULL;
	size_t ocalloc_size = sizeof(ms_sgx_oc_cpuidex_t);
	void *__tmp = NULL;

	void *__tmp_cpuinfo = NULL;

	CHECK_ENCLAVE_POINTER(cpuinfo, _len_cpuinfo);

	if (ADD_ASSIGN_OVERFLOW(ocalloc_size, (cpuinfo != NULL) ? _len_cpuinfo : 0))
		return SGX_ERROR_INVALID_PARAMETER;

	__tmp = sgx_ocalloc(ocalloc_size);
	if (__tmp == NULL) {
		sgx_ocfree();
		return SGX_ERROR_UNEXPECTED;
	}
	ms = (ms_sgx_oc_cpuidex_t*)__tmp;
	__tmp = (void *)((size_t)__tmp + sizeof(ms_sgx_oc_cpuidex_t));
	ocalloc_size -= sizeof(ms_sgx_oc_cpuidex_t);

	if (cpuinfo != NULL) {
		ms->ms_cpuinfo = (int*)__tmp;
		__tmp_cpuinfo = __tmp;
		if (_len_cpuinfo % sizeof(*cpuinfo) != 0) {
			sgx_ocfree();
			return SGX_ERROR_INVALID_PARAMETER;
		}
		memset(__tmp_cpuinfo, 0, _len_cpuinfo);
		__tmp = (void *)((size_t)__tmp + _len_cpuinfo);
		ocalloc_size -= _len_cpuinfo;
	} else {
		ms->ms_cpuinfo = NULL;
	}
	
	ms->ms_leaf = leaf;
	ms->ms_subleaf = subleaf;
	status = sgx_ocall(1, ms);

	if (status == SGX_SUCCESS) {
		if (cpuinfo) {
			if (memcpy_s((void*)cpuinfo, _len_cpuinfo, __tmp_cpuinfo, _len_cpuinfo)) {
				sgx_ocfree();
				return SGX_ERROR_UNEXPECTED;
			}
		}
	}
	sgx_ocfree();
	return status;
}

sgx_status_t SGX_CDECL sgx_thread_wait_untrusted_event_ocall(int* retval, const void* self)
{
	sgx_status_t status = SGX_SUCCESS;

	ms_sgx_thread_wait_untrusted_event_ocall_t* ms = NULL;
	size_t ocalloc_size = sizeof(ms_sgx_thread_wait_untrusted_event_ocall_t);
	void *__tmp = NULL;


	__tmp = sgx_ocalloc(ocalloc_size);
	if (__tmp == NULL) {
		sgx_ocfree();
		return SGX_ERROR_UNEXPECTED;
	}
	ms = (ms_sgx_thread_wait_untrusted_event_ocall_t*)__tmp;
	__tmp = (void *)((size_t)__tmp + sizeof(ms_sgx_thread_wait_untrusted_event_ocall_t));
	ocalloc_size -= sizeof(ms_sgx_thread_wait_untrusted_event_ocall_t);

	ms->ms_self = self;
	status = sgx_ocall(2, ms);

	if (status == SGX_SUCCESS) {
		if (retval) *retval = ms->ms_retval;
	}
	sgx_ocfree();
	return status;
}

sgx_status_t SGX_CDECL sgx_thread_set_untrusted_event_ocall(int* retval, const void* waiter)
{
	sgx_status_t status = SGX_SUCCESS;

	ms_sgx_thread_set_untrusted_event_ocall_t* ms = NULL;
	size_t ocalloc_size = sizeof(ms_sgx_thread_set_untrusted_event_ocall_t);
	void *__tmp = NULL;


	__tmp = sgx_ocalloc(ocalloc_size);
	if (__tmp == NULL) {
		sgx_ocfree();
		return SGX_ERROR_UNEXPECTED;
	}
	ms = (ms_sgx_thread_set_untrusted_event_ocall_t*)__tmp;
	__tmp = (void *)((size_t)__tmp + sizeof(ms_sgx_thread_set_untrusted_event_ocall_t));
	ocalloc_size -= sizeof(ms_sgx_thread_set_untrusted_event_ocall_t);

	ms->ms_waiter = waiter;
	status = sgx_ocall(3, ms);

	if (status == SGX_SUCCESS) {
		if (retval) *retval = ms->ms_retval;
	}
	sgx_ocfree();
	return status;
}

sgx_status_t SGX_CDECL sgx_thread_setwait_untrusted_events_ocall(int* retval, const void* waiter, const void* self)
{
	sgx_status_t status = SGX_SUCCESS;

	ms_sgx_thread_setwait_untrusted_events_ocall_t* ms = NULL;
	size_t ocalloc_size = sizeof(ms_sgx_thread_setwait_untrusted_events_ocall_t);
	void *__tmp = NULL;


	__tmp = sgx_ocalloc(ocalloc_size);
	if (__tmp == NULL) {
		sgx_ocfree();
		return SGX_ERROR_UNEXPECTED;
	}
	ms = (ms_sgx_thread_setwait_untrusted_events_ocall_t*)__tmp;
	__tmp = (void *)((size_t)__tmp + sizeof(ms_sgx_thread_setwait_untrusted_events_ocall_t));
	ocalloc_size -= sizeof(ms_sgx_thread_setwait_untrusted_events_ocall_t);

	ms->ms_waiter = waiter;
	ms->ms_self = self;
	status = sgx_ocall(4, ms);

	if (status == SGX_SUCCESS) {
		if (retval) *retval = ms->ms_retval;
	}
	sgx_ocfree();
	return status;
}

sgx_status_t SGX_CDECL sgx_thread_set_multiple_untrusted_events_ocall(int* retval, const void** waiters, size_t total)
{
	sgx_status_t status = SGX_SUCCESS;
	size_t _len_waiters = total * sizeof(void*);

	ms_sgx_thread_set_multiple_untrusted_events_ocall_t* ms = NULL;
	size_t ocalloc_size = sizeof(ms_sgx_thread_set_multiple_untrusted_events_ocall_t);
	void *__tmp = NULL;


	CHECK_ENCLAVE_POINTER(waiters, _len_waiters);

	if (ADD_ASSIGN_OVERFLOW(ocalloc_size, (waiters != NULL) ? _len_waiters : 0))
		return SGX_ERROR_INVALID_PARAMETER;

	__tmp = sgx_ocalloc(ocalloc_size);
	if (__tmp == NULL) {
		sgx_ocfree();
		return SGX_ERROR_UNEXPECTED;
	}
	ms = (ms_sgx_thread_set_multiple_untrusted_events_ocall_t*)__tmp;
	__tmp = (void *)((size_t)__tmp + sizeof(ms_sgx_thread_set_multiple_untrusted_events_ocall_t));
	ocalloc_size -= sizeof(ms_sgx_thread_set_multiple_untrusted_events_ocall_t);

	if (waiters != NULL) {
		ms->ms_waiters = (const void**)__tmp;
		if (_len_waiters % sizeof(*waiters) != 0) {
			sgx_ocfree();
			return SGX_ERROR_INVALID_PARAMETER;
		}
		if (memcpy_s(__tmp, ocalloc_size, waiters, _len_waiters)) {
			sgx_ocfree();
			return SGX_ERROR_UNEXPECTED;
		}
		__tmp = (void *)((size_t)__tmp + _len_waiters);
		ocalloc_size -= _len_waiters;
	} else {
		ms->ms_waiters = NULL;
	}
	
	ms->ms_total = total;
	status = sgx_ocall(5, ms);

	if (status == SGX_SUCCESS) {
		if (retval) *retval = ms->ms_retval;
	}
	sgx_ocfree();
	return status;
}

