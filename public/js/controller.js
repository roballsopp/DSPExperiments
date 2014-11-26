function addFile(evt) {
	
	var files = evt.target.files; // FileList object

	for ( var f = 0; f < files.length; f++ ) {
		console.log("Sending...");
		AsyncExec(decode, files[f], utils.isLittleEndian)
		.then(function (signal) {
			console.log(signal);
		}, function (error) {
			console.error();
		});
	}
	
}

realUploadFileBtn.addEventListener('change', addFile);

