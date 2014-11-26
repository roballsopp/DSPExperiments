var realUploadFileBtn = utils.createElement("input", {type: "file"});
var uploadFileBtn = document.getElementById("upload-file-btn");

uploadFileBtn.addEventListener('click', function (e) {
	realUploadFileBtn.dispatchEvent(new Event('click'));
});